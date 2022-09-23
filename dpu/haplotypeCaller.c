#include "haplotypeCaller.h"
#include "fixedComputation.h"
#include "./../host/constants.h"
#include "mutex.h"
#include <limits.h>


BARRIER_INIT(my_barrier, NR_TASKLETS);


//results 
__host uint64_t nb_cycles;
__mram_noinit int64_t likelihoods[MAX_HAPLOTYPE_NUM*MAX_READ_NUM];

//WRAM Data
__host uint64_t nr_reads;
__host uint64_t nr_haplotypes;
__host uint64_t dpu_inactive;
__host uint32_t haplotype_region_starts[MAX_REGIONS_PER_DPU+1];
__host uint32_t read_region_starts[MAX_REGIONS_PER_DPU+1];
  
//MRAM 
__mram_noinit uint64_t mram_reads_len[MAX_READ_NUM];
__mram_noinit char mram_reads_array[MAX_READ_NUM * MAX_READ_LENGTH];//[MAX_READ_NUM][MAX_READ_LENGTH];
__mram_noinit uint32_t mram_priors[2 * MAX_READ_NUM * MAX_READ_LENGTH];//[MAX_READ_NUM][2 * MAX_READ_LENGTH];
__mram_noinit uint64_t mram_haplotypes_len[MAX_HAPLOTYPE_NUM];
__mram_noinit int64_t mram_haplotypes_val[MAX_HAPLOTYPE_NUM];
__mram_noinit char mram_haplotypes_array[MAX_HAPLOTYPE_NUM * MAX_HAPLOTYPE_LENGTH];//[MAX_HAPLOTYPE_NUM][MAX_HAPLOTYPE_LENGTH];
__mram_noinit int32_t mram_matchToIndelArray[MAX_READ_NUM * MAX_READ_LENGTH];



//WRAM Data
//these transition values are correct for 9 decimal bits (otherwise it has to be recalculated)
int transition[TRANS_PROB_ARRAY_LENGTH] = { 0, -23, -2048, -512, -2048, -512, -2304 };
//fixed values corresponding to transition doubles { 0.999936759471893310546875, 0.999968377, 0.999968377, 0.9, 0.1, 0.1 };
__dma_aligned int64_t res[NR_TASKLETS];



__dma_aligned uint64_t reads_len[NR_TASKLETS];
__dma_aligned char reads_array[NR_TASKLETS][MAX_READ_LENGTH];
__dma_aligned int priors[NR_TASKLETS][2 * MAX_READ_LENGTH];

__dma_aligned uint64_t haplotypes_len_buffer[NR_WRAM_HAPLOTYPES]; //instead of sending haplotypes lengths now we will send the log10(1/hapLength) fixed number
__dma_aligned int64_t haplotypes_val_buffer[NR_WRAM_HAPLOTYPES];
uint32_t haplotypes_buffer_start;
uint32_t haplotypes_buffer_end;
__dma_aligned char haplotypes_buffer[NR_WRAM_HAPLOTYPES * MAX_HAPLOTYPE_LENGTH];



__dma_aligned int32_t matchToIndelArray[NR_TASKLETS][MAX_READ_LENGTH];


//Matrix cache (WRAM)
__dma_aligned int MATCH_CACHE[NR_TASKLETS][MATRIX_LINES][MAX_READ_LENGTH + 1];
__dma_aligned int INSERTION_CACHE[NR_TASKLETS][MAX_READ_LENGTH + 1];
__dma_aligned int DELETION_CACHE[NR_TASKLETS][MAX_READ_LENGTH + 1];

int fixedAdd(int a, int b) {
	if (a == INT_MIN || b == INT_MIN) {
		return INT_MIN;
	}
	int a1 = a;
	int b1 = b;

	int sum = a1 + b1;

	if (((~(a1 ^ b1) & (a1 ^ sum)) & INT_MIN) != 0) {
		sum = (a1 > 0 ? INT_MAX : INT_MIN);
	}
	return sum;
}

void initialize_matrices(uint32_t id, uint32_t haplotype_buffer_idx) {
	for (uint32_t idx = 0; idx < reads_len[id]; idx++) {
		MATCH_CACHE[id][0][idx] = INT_MIN;
		INSERTION_CACHE[id][idx] = INT_MIN;
		DELETION_CACHE[id][idx] = INT_MIN;
	}
	for (uint32_t idx = 0; idx < MATRIX_LINES; idx++) {
		MATCH_CACHE[id][idx][0] = INT_MIN;
		DELETION_CACHE[id][0] = haplotypes_val_buffer[haplotype_buffer_idx];
	}
}

void allocate_read_for_tasklet(uint32_t read_idx_th0, uint32_t tasklet_id) {
	mram_read(&mram_reads_len[read_idx_th0], &reads_len[tasklet_id], sizeof(uint64_t));
	int read_offset = read_idx_th0 * MAX_READ_LENGTH;
	int read_size = ((reads_len[tasklet_id] % 8 == 0) ? reads_len[tasklet_id] : reads_len[tasklet_id] + 8 - reads_len[tasklet_id] % 8) * sizeof(char);
	assert(read_size <= MAX_READ_LENGTH);
	mram_read(&mram_reads_array[read_offset], &reads_array[tasklet_id][0], read_size * sizeof(char));
	assert(2 * read_size * sizeof(int) <= 2 * MAX_READ_LENGTH * sizeof(uint32_t));
	mram_read(&mram_priors[2 * read_offset], &priors[tasklet_id][0], 2 * read_size * sizeof(int));
	mram_read(&mram_matchToIndelArray[read_offset], &matchToIndelArray[tasklet_id][0], read_size * sizeof(int32_t));
}
 

void allocate_haplotypes() {
	//IMPORTANT: the transfer size must be a multiple of 8, we need to make a correction in some cases
	//We need to keep MAX_HAPLOTYPE_LENGTH multiple of 4 
	int transfer_size = nr_haplotypes + (nr_haplotypes % 2 == 1);
	assert(transfer_size <= MAX_HAPLOTYPE_NUM);
	mram_read(mram_haplotypes_len, haplotypes_len_buffer, transfer_size * sizeof(uint64_t));
	mram_read(mram_haplotypes_val, haplotypes_val_buffer, transfer_size * sizeof(int64_t));
	int start;
	for (start=0; start<transfer_size*MAX_HAPLOTYPE_LENGTH*sizeof(char) - LIMIT; start+=LIMIT) {
		mram_read((__mram_ptr void*)mram_haplotypes_array + start, (void*)haplotypes_buffer + start, LIMIT);
	}
	mram_read((__mram_ptr void*)mram_haplotypes_array + start, (void*)haplotypes_buffer + start, transfer_size * MAX_HAPLOTYPE_LENGTH * sizeof(char) - start);
}


uint32_t free_read_idx;
int32_t last_region_allocated;
uint32_t current_region[NR_TASKLETS];
MUTEX_INIT(task_reservation_mutex);

uint32_t reserve_read(int tasklet_id) {
	current_region[tasklet_id] = MAX_REGIONS_PER_DPU+1;

	mutex_lock(task_reservation_mutex);

	uint32_t result = free_read_idx++;
	// If the new read is in a new region, allocate necessary haplotypes
	if (result < nr_reads && result>=read_region_starts[last_region_allocated+1]) {
		last_region_allocated++;
		uint32_t number_of_haplotypes = haplotype_region_starts[last_region_allocated+1]-haplotype_region_starts[last_region_allocated];
		assert(number_of_haplotypes < NR_WRAM_HAPLOTYPES);
		// If there isn't enough room in the haplotypes buffer,
		// Actively wait for other tasklets to finish
		// Active wait isn't optimal but hopefully it shouldn't happen too often.
		// TODO: investigate cost of active wait.
		while (haplotypes_buffer_start != haplotypes_buffer_end && number_of_haplotypes > (haplotypes_buffer_start-haplotypes_buffer_end)%NR_WRAM_HAPLOTYPES) {
			int min_region = last_region_allocated;
			for (int id=0; id<NR_TASKLETS; id++) {
				if (current_region[id] < min_region) {
					min_region = current_region[id];
				}
			}
			// update the haplotypes_buffer_start in case any thread finished processing a region.
			assert(min_region<NR_REGIONS);
			haplotypes_buffer_start = haplotype_region_starts[min_region]%NR_WRAM_HAPLOTYPES;
		}

		int transfer_size;
		// Only copy as far as you can at first
		if (number_of_haplotypes+haplotypes_buffer_end > NR_WRAM_HAPLOTYPES) {
			assert(number_of_haplotypes<=NR_WRAM_HAPLOTYPES);
			transfer_size = (NR_WRAM_HAPLOTYPES-haplotypes_buffer_end);
		} else {
			transfer_size = number_of_haplotypes;
		}
		// Start copy at the end of the buffer.
		int haplotype_copy_start_index = haplotype_region_starts[last_region_allocated];
		mram_read(&mram_haplotypes_len[haplotype_copy_start_index], &haplotypes_len_buffer[haplotypes_buffer_end], transfer_size * sizeof(uint64_t));
		mram_read(&mram_haplotypes_val[haplotype_copy_start_index], &haplotypes_val_buffer[haplotypes_buffer_end], transfer_size * sizeof(uint64_t));
		// If everything couldn't be copied, go back to the start of the circular buffer and copy what's left
		if (number_of_haplotypes+haplotypes_buffer_start > NR_WRAM_HAPLOTYPES) {
			haplotype_copy_start_index += transfer_size;// Skip what was already copied
			transfer_size = (number_of_haplotypes-NR_WRAM_HAPLOTYPES+haplotypes_buffer_start);
			mram_read(&mram_haplotypes_len[haplotype_copy_start_index], haplotypes_len_buffer, transfer_size * sizeof(uint64_t));
			mram_read(&mram_haplotypes_val[haplotype_copy_start_index], haplotypes_val_buffer, transfer_size * sizeof(uint64_t));
		}
		// Now copy the haplotypes while making sure no transfer is bigger than LIMIT
		#define BUFFER_SIZE (MAX_HAPLOTYPE_LENGTH*NR_WRAM_HAPLOTYPES*sizeof(char))
		int transfer_end = (haplotypes_buffer_end+transfer_size)*MAX_HAPLOTYPE_LENGTH;
		for (int written = 0; haplotypes_buffer_end*MAX_HAPLOTYPE_LENGTH+written < transfer_end; ) {
			int write_length;
			if (written+LIMIT+haplotypes_buffer_end*MAX_HAPLOTYPE_LENGTH >= transfer_end) {
				write_length = transfer_end - (haplotypes_buffer_end*MAX_HAPLOTYPE_LENGTH*sizeof(char) + written);
			} else {
				write_length = LIMIT;
			}
			if ((haplotypes_buffer_end*MAX_HAPLOTYPE_LENGTH+written+write_length)%BUFFER_SIZE < haplotypes_buffer_end*MAX_HAPLOTYPE_LENGTH+written+write_length) {
				write_length = BUFFER_SIZE - (haplotypes_buffer_end*MAX_HAPLOTYPE_LENGTH*sizeof(char) + written);
			}
			mram_read((__mram_ptr void*)&mram_haplotypes_array[haplotype_copy_start_index*MAX_HAPLOTYPE_LENGTH*sizeof(char)] + written,
					  (void*)haplotypes_buffer + (haplotypes_buffer_end*MAX_HAPLOTYPE_LENGTH*sizeof(char) + written)%BUFFER_SIZE,
					  write_length);
			written += write_length;
		}
		// Update the position of buffer_end
		haplotypes_buffer_end += number_of_haplotypes;
		haplotypes_buffer_end %= NR_WRAM_HAPLOTYPES;
	}
	current_region[tasklet_id] = last_region_allocated;

	mutex_unlock(task_reservation_mutex);

	return result;
}


int main() {
	if (dpu_inactive == 1) { return 1; }
	// if (dpu_inactive != 1) { return 0; }// FIXME : remove that line
	thread_id_t tasklet_id = me();
	current_region[tasklet_id] = 0;
	// uint32_t rounds = nr_reads / NR_TASKLETS + (nr_reads % NR_TASKLETS != 0);
	if (tasklet_id == 0) {
		free_read_idx = 0;
		nb_cycles = 0;
		perfcounter_config(COUNT_CYCLES, true);
		last_region_allocated = -1;
		free_read_idx = 0;
		haplotypes_buffer_start = 0;
		haplotypes_buffer_end = 0;
	}
	/*
	if (tasklet_id == 0) {
		allocate_haplotypes();
	} 
	*/
	barrier_wait(&my_barrier);

	//for (uint32_t round = 0; round < rounds; round++) {
		//uint32_t	read_idx = tasklet_id + NR_TASKLETS * round;

	uint32_t read_idx;
	while ((read_idx = reserve_read(tasklet_id)) < nr_reads) {

		if (read_idx < nr_reads) {
			allocate_read_for_tasklet(read_idx, tasklet_id);
			uint32_t j;
			//From here each tasklet is in charge of reads tasklet_id, tasklet_id + NR_TASKLETS, tasklet_id + 2*NR_TASKLETS ...
			for (uint32_t haplotype_idx = haplotype_region_starts[current_region[tasklet_id]];
						  haplotype_idx < haplotype_region_starts[current_region[tasklet_id]+1];
						  haplotype_idx++) {
				uint32_t haplotype_buffer_idx = haplotype_idx%NR_WRAM_HAPLOTYPES;

				initialize_matrices(tasklet_id, haplotype_buffer_idx);
				res[tasklet_id] = INT_MIN;
				for (uint32_t i = 1; i <= haplotypes_len_buffer[haplotype_buffer_idx]; i++) {
					uint32_t indI = i % MATRIX_LINES;
					uint32_t indI0 = (i - 1) % MATRIX_LINES;
                    int insertion_diagonal_value = INSERTION_CACHE[tasklet_id][0];
                    int deletion_diagonal_value  = DELETION_CACHE[tasklet_id][0];
					for (j = 1; j <= reads_len[tasklet_id]-1; j++) {
						int prior;

						if (reads_array[tasklet_id][j - 1] == haplotypes_buffer[haplotype_buffer_idx * MAX_HAPLOTYPE_LENGTH + (i - 1)] ||
							reads_array[tasklet_id][j-1] == 'N' || haplotypes_buffer[haplotype_buffer_idx * MAX_HAPLOTYPE_LENGTH + (i - 1)] == 'N') {

							prior = priors[tasklet_id][2 * (j - 1)];
						} 
						else {
							prior = priors[tasklet_id][2 * (j - 1) + 1];
						}

						
						int matchToIndel = (j==1)? -2048:matchToIndelArray[tasklet_id][j-2];
						MATCH_CACHE[tasklet_id][indI][j] = fixedAdd(prior, log10SumLog10(log10SumLog10(fixedAdd(MATCH_CACHE[tasklet_id][indI0][j - 1], transition[matchToMatch]),
							fixedAdd(insertion_diagonal_value, transition[indelToMatch])),
							fixedAdd(deletion_diagonal_value, transition[indelToMatch])));

                        insertion_diagonal_value = INSERTION_CACHE[tasklet_id][j];
						INSERTION_CACHE[tasklet_id][j] = log10SumLog10(fixedAdd(MATCH_CACHE[tasklet_id][indI][j - 1], matchToIndel), fixedAdd(INSERTION_CACHE[tasklet_id][j - 1], transition[insertionToInsertion]));

                        deletion_diagonal_value = DELETION_CACHE[tasklet_id][j];
						DELETION_CACHE[tasklet_id][j] = log10SumLog10(fixedAdd(MATCH_CACHE[tasklet_id][indI0][j], matchToIndel), fixedAdd(DELETION_CACHE[tasklet_id][j], transition[deletionToDeletion]));
						
					}
					//The last iteration on the read length is quite different (different transition probability used in 
					//INSERTION and DELETION matrix element computation)
					//So we compute these last elements separately
					{
						int prior;

						if (reads_array[tasklet_id][j - 1] == haplotypes_buffer[haplotype_buffer_idx * MAX_HAPLOTYPE_LENGTH + (i - 1)] ||
							reads_array[tasklet_id][j-1] == 'N' || haplotypes_buffer[haplotype_buffer_idx * MAX_HAPLOTYPE_LENGTH + (i - 1)] == 'N') {

							prior = priors[tasklet_id][2 * (j - 1)];
						} 
						else {
							prior = priors[tasklet_id][2 * (j - 1) + 1];
						}

						
						int matchToIndel = (j==1)? -2048:matchToIndelArray[tasklet_id][j-2];
						MATCH_CACHE[tasklet_id][indI][j] = fixedAdd(prior, log10SumLog10(log10SumLog10(fixedAdd(MATCH_CACHE[tasklet_id][indI0][j - 1], transition[matchToMatch]),
							fixedAdd(insertion_diagonal_value, transition[indelToMatch])),
							fixedAdd(deletion_diagonal_value, transition[indelToMatch])));
						INSERTION_CACHE[tasklet_id][j] = log10SumLog10(fixedAdd(MATCH_CACHE[tasklet_id][indI][j - 1], transition[lastBaseTransition]), fixedAdd(INSERTION_CACHE[tasklet_id][j - 1], transition[insertionToInsertion]));
						DELETION_CACHE[tasklet_id][j] = log10SumLog10(fixedAdd(MATCH_CACHE[tasklet_id][indI0][j], transition[lastBaseTransition]), fixedAdd(DELETION_CACHE[tasklet_id][j], transition[deletionToDeletion]));
					}
					res[tasklet_id] = log10SumLog10(res[tasklet_id], log10SumLog10(MATCH_CACHE[tasklet_id][i % MATRIX_LINES][j/*reads_len[tasklet_id]*/], INSERTION_CACHE[tasklet_id][j/*reads_len[tasklet_id] */ ]));
					
				}

				mram_write(&res[tasklet_id], &likelihoods[haplotype_idx* MAX_READ_NUM + read_idx], sizeof(res[0]));
			}
		}

		
	}

	barrier_wait(&my_barrier);
	if (tasklet_id == 0) {
		nb_cycles = perfcounter_get();
	}
	barrier_wait(&my_barrier);
	return 1;
}
