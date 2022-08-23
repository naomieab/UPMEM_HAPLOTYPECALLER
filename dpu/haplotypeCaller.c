#include "haplotypeCaller.h"
#include "fixedComputation.h"
#include "./../host/constants.h"
#include "mutex.h"
#include <limits.h>
#include <stdio.h>


BARRIER_INIT(my_barrier, NR_TASKLETS);


//results 
__host uint64_t nb_cycles;
__mram_noinit int64_t likelihoods[MAX_HAPLOTYPE_NUM][MAX_READ_NUM];

//WRAM Data
__host uint64_t nr_reads;
__host uint64_t nr_haplotypes;
__host uint64_t dpu_inactive;
  
//MRAM 
__mram_noinit uint64_t mram_reads_len[MAX_READ_NUM];
__mram_noinit char mram_reads_array[MAX_READ_NUM * MAX_READ_LENGTH];//[MAX_READ_NUM][MAX_READ_LENGTH];
__mram_noinit uint32_t mram_priors[2 * MAX_READ_NUM * MAX_READ_LENGTH];//[MAX_READ_NUM][2 * MAX_READ_LENGTH];
__mram_noinit uint32_t mram_haplotypes_len[MAX_HAPLOTYPE_NUM];
__mram_noinit uint32_t mram_haplotypes_val[MAX_HAPLOTYPE_NUM];
__mram_noinit char mram_haplotypes_array[MAX_HAPLOTYPE_NUM * MAX_HAPLOTYPE_LENGTH];//[MAX_HAPLOTYPE_NUM][MAX_HAPLOTYPE_LENGTH];




//WRAM Data
//these transition values are correct for 9 decimal bits (otherwise it has to be recalculated)
int transition[TRANS_PROB_ARRAY_LENGTH] = { 0, -23, -2048, -512, -2048, -512, -2304 };
//fixed values corresponding to transition doubles { 0.999936759471893310546875, 0.999968377, 0.999968377, 0.9, 0.1, 0.1 };
__dma_aligned int64_t res[NR_TASKLETS];


__dma_aligned uint64_t reads_len[NR_TASKLETS];
__dma_aligned char reads_array[NR_TASKLETS][MAX_READ_LENGTH];
__dma_aligned int priors[NR_TASKLETS][2 * MAX_READ_LENGTH];
__dma_aligned uint32_t haplotypes_len[MAX_HAPLOTYPE_NUM]; //instead of sending haplotypes lengths now we will send the log10(1/hapLength) fixed number
__dma_aligned uint32_t haplotypes_val[MAX_HAPLOTYPE_NUM];
__dma_aligned char haplotypes_array[MAX_HAPLOTYPE_NUM][MAX_HAPLOTYPE_LENGTH];


//Matrix cache (WRAM)
__dma_aligned int MATCH_CACHE[NR_TASKLETS][MATRIX_LINES][MAX_READ_LENGTH + 1];
__dma_aligned int INSERTION_CACHE[NR_TASKLETS][MATRIX_LINES][MAX_READ_LENGTH + 1];
__dma_aligned int DELETION_CACHE[NR_TASKLETS][MATRIX_LINES][MAX_READ_LENGTH + 1];

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

void initialize_matrices(uint32_t id, uint32_t haplotype_idx) {
	for (uint32_t idx = 0; idx < reads_len[id]; idx++) {
		MATCH_CACHE[id][0][idx] = INT_MIN;
		INSERTION_CACHE[id][0][idx] = INT_MIN;
		DELETION_CACHE[id][0][idx] = INT_MIN;
	}
	for (uint32_t idx = 0; idx < MATRIX_LINES; idx++) {
		MATCH_CACHE[id][idx][0] = INT_MIN;
		INSERTION_CACHE[id][idx][0] = INT_MIN;
		DELETION_CACHE[id][idx][0] = haplotypes_val[haplotype_idx];
	}
}

void allocate_read_for_tasklet(uint32_t read_idx_th0, uint32_t tasklet_id) {
	mram_read(&mram_reads_len[read_idx_th0], &reads_len[tasklet_id], sizeof(uint64_t));
	int read_offset = read_idx_th0 * MAX_READ_LENGTH;
	int read_size = ((reads_len[tasklet_id] % 8 == 0) ? reads_len[tasklet_id] : reads_len[tasklet_id] + 8 - reads_len[tasklet_id] % 8) * sizeof(char);
	assert(read_size <= MAX_READ_LENGTH);
	mram_read(&mram_reads_array[read_offset], &reads_array[tasklet_id][0], read_size * sizeof(char));
	assert(2*read_size*sizeof(int) <= 2*MAX_READ_LENGTH*sizeof(uint32_t));
	mram_read(&mram_priors[2*read_offset], &priors[tasklet_id][0], 2 * read_size * sizeof(int));
}


void allocate_haplotypes() {
	//IMPORTANT: the transfer size must be a multiple of 8, we need to make a correction in some cases
	//We need to keep MAX_HAPLOTYPE_LENGTH multiple of 4 
	int transfer_size = nr_haplotypes + (nr_haplotypes % 2 == 1);
	assert(transfer_size <= MAX_HAPLOTYPE_NUM);
	mram_read(mram_haplotypes_len, haplotypes_len, transfer_size * sizeof(uint32_t));
	mram_read(mram_haplotypes_val, haplotypes_val, transfer_size * sizeof(uint32_t));
	mram_read(mram_haplotypes_array, haplotypes_array, transfer_size * MAX_HAPLOTYPE_LENGTH * sizeof(char));
}


uint32_t free_read_idx;
MUTEX_INIT(task_reservation_mutex);
uint32_t reserve_read(int tasklet_id) {
	mutex_lock(task_reservation_mutex);
	uint32_t result = free_read_idx++;
	mutex_unlock(task_reservation_mutex);
	return result;
}



int main() {
	if (dpu_inactive == 1) { return 1; }
	thread_id_t tasklet_id = me();
	uint32_t rounds = nr_reads / NR_TASKLETS + (nr_reads % NR_TASKLETS != 0);
	if (tasklet_id == 0) {
		free_read_idx = 0;
		nb_cycles = 0;
		perfcounter_config(COUNT_CYCLES, true);
	}
	if (tasklet_id == 0) {
		allocate_haplotypes();
	} 
	barrier_wait(&my_barrier);
	
	//for (uint32_t round = 0; round < rounds; round++) {
		//uint32_t	read_idx = tasklet_id + NR_TASKLETS * round;

	uint32_t read_idx;
	while ((read_idx = reserve_read(tasklet_id)) < nr_reads) {
		
		if (read_idx < nr_reads) {
			allocate_read_for_tasklet(read_idx, tasklet_id);
			uint32_t j;
			//From here each tasklet is in charge of reads tasklet_id, tasklet_id + NR_TASKLETS, tasklet_id + 2*NR_TASKLETS ...
			for (uint32_t haplotype_idx = 0; haplotype_idx < nr_haplotypes; haplotype_idx++) {

				initialize_matrices(tasklet_id, haplotype_idx);
				res[tasklet_id] = INT_MIN;
				for (uint32_t i = 1; i <= haplotypes_len[haplotype_idx]; i++) {
					uint32_t indI = i % MATRIX_LINES;
					uint32_t indI0 = (i - 1) % MATRIX_LINES;
					for (j = 1; j <= reads_len[tasklet_id]; j++) {
						int prior;
						if (reads_array[tasklet_id][j - 1] == haplotypes_array[haplotype_idx][i - 1]) {
							prior = priors[tasklet_id][2 * (j - 1)];
						}
						else {
							prior = priors[tasklet_id][2 * (j - 1) + 1];
						}

						MATCH_CACHE[tasklet_id][indI][j] = fixedAdd(prior, log10SumLog10(log10SumLog10(fixedAdd(MATCH_CACHE[tasklet_id][indI0][j - 1], transition[matchToMatch]),
							fixedAdd(INSERTION_CACHE[tasklet_id][indI0][j - 1], transition[indelToMatch])),
							fixedAdd(DELETION_CACHE[tasklet_id][indI0][j - 1], transition[indelToMatch])));

						INSERTION_CACHE[tasklet_id][indI][j] = log10SumLog10(fixedAdd(MATCH_CACHE[tasklet_id][indI][j - 1], transition[matchToInsertion]), fixedAdd(INSERTION_CACHE[tasklet_id][indI][j - 1], transition[insertionToInsertion]));

						DELETION_CACHE[tasklet_id][indI][j] = log10SumLog10(fixedAdd(MATCH_CACHE[tasklet_id][indI0][j], transition[matchToDeletion]), fixedAdd(DELETION_CACHE[tasklet_id][indI0][j], transition[deletionToDeletion]));

					}
					//The last iteration on the read length is quite different (different transition probability used in 
					//INSERTION and DELETION matrix element computation)
					//So we recompute these two last elements
					j--; //recompute at the last index
					INSERTION_CACHE[tasklet_id][indI][j] = log10SumLog10(fixedAdd(MATCH_CACHE[tasklet_id][indI][j - 1], transition[lastBaseTransition]), fixedAdd(INSERTION_CACHE[tasklet_id][indI][j - 1], transition[insertionToInsertion]));
					DELETION_CACHE[tasklet_id][indI][j] = log10SumLog10(fixedAdd(MATCH_CACHE[tasklet_id][indI0][j], transition[lastBaseTransition]), fixedAdd(DELETION_CACHE[tasklet_id][indI0][j], transition[deletionToDeletion]));
					res[tasklet_id] = log10SumLog10(res[tasklet_id], log10SumLog10(MATCH_CACHE[tasklet_id][i % MATRIX_LINES][reads_len[tasklet_id] - 1], INSERTION_CACHE[tasklet_id][i % MATRIX_LINES][reads_len[tasklet_id] - 1]));
				}

				mram_write(&res[tasklet_id], &likelihoods[haplotype_idx][read_idx], sizeof(res[0]));
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
