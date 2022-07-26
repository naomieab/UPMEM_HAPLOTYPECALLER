#include "haplotypeCaller.h"
#include "fixedComputation.h"
#include "./../host/constants.h"
#include <limits.h>


BARRIER_INIT(my_barrier, NR_TASKLETS);

//results
__host uint32_t nb_cycles;
__mram_noinit int likelihoods[MAX_HAPLOTYPE_NUM][MAX_READ_NUM];


//WRAM Data
__host uint32_t nr_reads;
__host uint32_t nr_haplotypes;

//MRAM
__mram_noinit uint32_t mram_reads_len[MAX_READ_NUM];
__mram_noinit char mram_reads_array[MAX_READ_NUM][MAX_READ_LENGTH];
__mram_noinit uint32_t mram_priors[MAX_READ_NUM][2 * MAX_READ_LENGTH];
__mram_noinit uint32_t mram_haplotypes_len[MAX_HAPLOTYPE_NUM];
__mram_noinit uint32_t mram_haplotypes_val[MAX_HAPLOTYPE_NUM];
__mram_noinit char mram_haplotypes_array[MAX_HAPLOTYPE_NUM][MAX_HAPLOTYPE_LENGTH];




//WRAM Data
int transition[TRANS_PROB_ARRAY_LENGTH] = { 511, 511, 511, 460, 51, 51 }; //fixed values corresponding to transition doubles { 0.999936759471893310546875, 0.999968377, 0.999968377, 0.9, 0.1, 0.1 };
int res[NR_TASKLETS];


uint32_t reads_len[NR_TASKLETS];
char reads_array[NR_TASKLETS][MAX_READ_LENGTH];
int priors[NR_TASKLETS][2 * MAX_READ_LENGTH];
uint32_t haplotypes_len[MAX_HAPLOTYPE_NUM]; //instead of sending haplotypes lengths now we will send the log10(1/hapLength) fixed number
uint32_t haplotypes_val[MAX_HAPLOTYPE_NUM];
char haplotypes_array[MAX_HAPLOTYPE_NUM][MAX_HAPLOTYPE_LENGTH];


//Matrix cache (WRAM)
int MATCH_CACHE[NR_TASKLETS][MATRIX_LINES][MAX_READ_LENGTH + 1];
int INSERTION_CACHE[NR_TASKLETS][MATRIX_LINES][MAX_READ_LENGTH + 1];
int DELETION_CACHE[NR_TASKLETS][MATRIX_LINES][MAX_READ_LENGTH + 1];

void initialize_matrices(uint32_t id, uint32_t haplotype_idx) {
	for (uint32_t idx = 0; idx < reads_len[id]; idx++) {
		MATCH_CACHE[id][0][idx] = INT_MIN;
		INSERTION_CACHE[id][0][idx] = INT_MIN;
		DELETION_CACHE[id][0][idx] = INT_MIN;
	}

	for (uint32_t idx = 1; idx < MATRIX_LINES; idx++) {
		MATCH_CACHE[id][idx][0] = INT_MIN;
		INSERTION_CACHE[id][idx][0] = haplotypes_val[haplotype_idx];
		DELETION_CACHE[id][idx][0] = INT_MIN;
	}
}

void allocate_reads(uint32_t read_idx_th0, uint32_t active_threads) {

	mram_read(&mram_reads_len[read_idx_th0], reads_len, NR_TASKLETS * sizeof(uint32_t));
	for (int i = 0; i < active_threads; i++) {
		int read_size = ((reads_len[i] % 8 == 0) ? reads_len[i] : reads_len[i] + 8 - reads_len[i] % 8) * sizeof(char);
		mram_read(&mram_reads_array[read_idx_th0 + i][0], &reads_array[i][0], read_size * sizeof(char));
		mram_read(&mram_priors[read_idx_th0 + i][0], &priors[i][0], 2 * read_size * sizeof(int));
	}
}


void allocate_haplotypes() {
	//IMPORTANT: the transfer size must be a multiple of 8, we need to make a correction in some cases
	//We need to keep MAX_HAPLOTYPE_LENGTH multiple of 4 
	int transfer_size = nr_haplotypes + (nr_haplotypes % 2 == 1);
	mram_read(mram_haplotypes_len, haplotypes_len, transfer_size * sizeof(uint32_t));
	mram_read(mram_haplotypes_val, haplotypes_val, transfer_size * sizeof(uint32_t));
	mram_read(mram_haplotypes_array, haplotypes_array, transfer_size * MAX_HAPLOTYPE_LENGTH * sizeof(char));

}



int main() {
	thread_id_t tasklet_id = me();
	uint32_t rounds = nr_reads / NR_TASKLETS + (nr_reads % NR_TASKLETS != 0);
	if (tasklet_id == 0) {
		nb_cycles = 0;
		perfcounter_config(COUNT_CYCLES, true);
	}
	if (tasklet_id == 0) {
		allocate_haplotypes();
	}
	barrier_wait(&my_barrier);

	for (uint32_t round = 0; round < rounds; round++) {
		uint32_t	read_idx = tasklet_id + NR_TASKLETS * round;
		
		if (tasklet_id == 0) {
			if (round == rounds - 1) {
				allocate_reads(read_idx, nr_reads - read_idx);
			}
			else {
				allocate_reads(read_idx, NR_TASKLETS);
			}
		}
		barrier_wait(&my_barrier);

		//From here each tasklet is on charge of reads tasklet_id, tasklet_id + NR_TASKLET, tasklet_id + 2*NR_TASKLET ...
		for (uint32_t haplotype_idx = 0; haplotype_idx < nr_haplotypes; haplotype_idx++) {

			if (read_idx < nr_reads) {

				initialize_matrices(tasklet_id, haplotype_idx);
				res[tasklet_id] = 0;

				for (uint32_t i = 1; i <= haplotypes_len[haplotype_idx]; i++) {

					for (uint32_t j = 1; j <= reads_len[tasklet_id]; j++) {

						int prior;
						if (reads_array[tasklet_id][j - 1] == haplotypes_array[haplotype_idx][i - 1]) {
							prior = priors[tasklet_id][2 * (j - 1)];
						}
						else {
							prior = priors[tasklet_id][2 * (j - 1) + 1];
						}

						uint32_t indI = i % MATRIX_LINES;
						uint32_t indI0 = (i - 1) % MATRIX_LINES;

						MATCH_CACHE[tasklet_id][indI][j] = fixedAdd(prior, log10SumLog10(log10SumLog10(fixedAdd(MATCH_CACHE[tasklet_id][indI0][j - 1], transition[matchToMatch]),
							fixedAdd(INSERTION_CACHE[tasklet_id][indI0][j - 1], transition[indelToMatch])),
							fixedAdd(DELETION_CACHE[tasklet_id][indI0][j - 1], transition[indelToMatch])));



						INSERTION_CACHE[tasklet_id][indI][j] = log10SumLog10(fixedAdd(MATCH_CACHE[tasklet_id][indI][j - 1], transition[matchToInsertion]), fixedAdd(INSERTION_CACHE[tasklet_id][indI][j - 1], transition[insertionToInsertion]));





						DELETION_CACHE[tasklet_id][indI][j] = log10SumLog10(fixedAdd(MATCH_CACHE[tasklet_id][indI0][j], transition[matchToDeletion]), fixedAdd(DELETION_CACHE[tasklet_id][indI0][j], transition[deletionToDeletion]));


					}

					res[tasklet_id] = log10SumLog10(res[tasklet_id], log10SumLog10(MATCH_CACHE[tasklet_id][i % MATRIX_LINES][reads_len[tasklet_id] - 1], INSERTION_CACHE[tasklet_id][i % MATRIX_LINES][reads_len[tasklet_id] - 1]));

				}

			}

			barrier_wait(&my_barrier);

			if (tasklet_id == 0) {
				mram_write(&res, &likelihoods[haplotype_idx][read_idx], sizeof(res));
			}

			barrier_wait(&my_barrier);
		}


		barrier_wait(&my_barrier);
	}

	if (tasklet_id == 0) {
		nb_cycles = perfcounter_get();
	}
	barrier_wait(&my_barrier);
	return 1;
}
