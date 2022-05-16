#include "haplotypeCaller.h"
#include "fixedComputation.h"
#include <limits.h>

BARRIER_INIT(my_barrier, NR_TASKLETS);

__host uint32_t nr_reads = 1;
__host uint32_t nr_haplotypes = 2;
/**
* Order in the heap is defined as :
* 0. reads_len,
* 1. reads_arr
* 2. prior_arr
* 3. hap_len
* 4. hap_arr
* 5. hap_val
*/
__host uint32_t heap_offsets[6];

__host int result[8];
__mram_noinit int likelihoods[MAX_HAPLOTYPES][MAX_READS];


//__host double read_check[50];
//__host double prior_check[238];
//__host double read_check2[120];
//__host double prior_check2[238];


//MRAM variables
//all the variables are passed through the heap and the result array is declared in MRAM


int transition[TRANS_PROB_ARRAY_LENGTH] = { 511, 511, 511, 460, 51, 51 }; //fixed values corresponding to transition doubles { 0.999936759471893310546875, 0.999968377, 0.999968377, 0.9, 0.1, 0.1 };
int res[NR_TASKLETS];


//memory in the WRAM, will be allocated with mem alloc and accessed by all threads
//We keep in WRAM 16 reads and their qualities, and all the haplotypes
uint32_t reads_len[NR_TASKLETS];
char* reads_arr[NR_TASKLETS];
int* prior_arr[NR_TASKLETS];
uint32_t* haplotypes_len; //instead of sending haplotypes lengths now we will send the log10(1/hapLength) fixed number
uint32_t* haplotypes_val;
char** haplotypes_arr;



int* MATCH_CACHE[NR_TASKLETS][MATRIX_LINES];
int* INSERTION_CACHE[NR_TASKLETS][MATRIX_LINES];
int* DELETION_CACHE[NR_TASKLETS][MATRIX_LINES];

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


/**
* n_reads is the number of reads to allocate (NR_TASKLET in all rounds except the last one that can require less)
*/
void allocate_reads(uint32_t round) {
	uintptr_t heap_address = ((uintptr_t)DPU_MRAM_HEAP_POINTER);
	uint32_t active_threads = NR_TASKLETS;//(round == nr_reads / NR_TASKLET) ? nr_reads%NR_TASKLET: NR_TASKLET;
	uint32_t i = 0;

	mram_read((__mram_ptr void*)(heap_address + heap_offsets[0]), reads_len, sizeof(reads_len));
	heap_offsets[0] += sizeof(reads_len);

	
	for (i = 0; i < active_threads; i++) {
		int read_size = ((reads_len[i] % 8 == 0) ? reads_len[i] : reads_len[i] + 8 - reads_len[i] % 8) * sizeof(char);
		reads_arr[i] = mem_alloc(read_size);
		mram_read((__mram_ptr void*)(heap_address + heap_offsets[2]), reads_arr[i], read_size);
		heap_offsets[2] += read_size;
	}
	for (i = 0; i < active_threads; i++) {
		int prior_size = 2 * (reads_len[i]) * sizeof(int);
		prior_arr[i] = mem_alloc(prior_size);
		mram_read((__mram_ptr void*)(heap_address + heap_offsets[3]), prior_arr[i], prior_size);
		heap_offsets[3] += prior_size;
	}
}

void allocate_haplotypes() {
	uintptr_t heap_address = ((uintptr_t)DPU_MRAM_HEAP_POINTER);
	uint32_t haplotypes_len_size = (nr_haplotypes + nr_haplotypes % 2) * sizeof(uint32_t);
	uint32_t tmp = heap_offsets[4];

	haplotypes_len = mem_alloc(haplotypes_len_size);
	mram_read((__mram_ptr void*)(heap_address + heap_offsets[1]), haplotypes_len, haplotypes_len_size);



	
	haplotypes_arr = mem_alloc(nr_haplotypes * sizeof(char*));

	for (uint32_t i = 0; i < nr_haplotypes; i++) {
		int hap_size = ((haplotypes_len[i] % 8 == 0) ? haplotypes_len[i] : haplotypes_len[i] + 8 - haplotypes_len[i] % 8) * sizeof(char);
		haplotypes_arr[i] = mem_alloc(hap_size);
		mram_read((__mram_ptr void*)(heap_address + heap_offsets[4]), haplotypes_arr[i], hap_size);
		heap_offsets[4] += hap_size;
	}
	heap_offsets[4] = tmp; //(because we want to read same haplotypes every time)
	
	haplotypes_val = mem_alloc(haplotypes_len_size);
	mram_read((__mram_ptr void*)(heap_address + heap_offsets[5]), haplotypes_val, haplotypes_len_size);

}


int main() {
	thread_id_t tasklet_id = me();
	uint32_t rounds = nr_reads / NR_TASKLETS + (nr_reads % NR_TASKLETS != 0);

	for (uint32_t round = 0; round < rounds; round++) {

		if (tasklet_id == 0) {
			mem_reset();
			allocate_haplotypes();
			allocate_reads(round);
		}

		uint32_t	read_idx = tasklet_id + NR_TASKLETS * round;
		barrier_wait(&my_barrier);

		if (tasklet_id == 0) {
			//result[1] = 1;
		}
		//From here each tasklet is on charge of reads tasklet_id, tasklet_id + NR_TASKLET, tasklet_id + 2*NR_TASKLET ...

		if (read_idx <= nr_reads) {
			//First allocate matrix according to current read length
			for (uint32_t i = 0; i < MATRIX_LINES; i++) {
				MATCH_CACHE[tasklet_id][i] = mem_alloc((reads_len[tasklet_id] + 1) * sizeof(int));
				INSERTION_CACHE[tasklet_id][i] = mem_alloc((reads_len[tasklet_id] + 1) * sizeof(int));
				DELETION_CACHE[tasklet_id][i] = mem_alloc((reads_len[tasklet_id] + 1) * sizeof(int));
			}
			for (uint32_t haplotype_idx = 0; haplotype_idx < nr_haplotypes; haplotype_idx++) {
				if (tasklet_id == 0) {
					//result[2] = 1;
				}
				initialize_matrices(tasklet_id, haplotype_idx);

				res[tasklet_id] = 0;
				if (tasklet_id == 0) {
					result[3] = 1;
				}
				
				for (uint32_t i = 1; i <= haplotypes_len[haplotype_idx]; i++) {
					for (uint32_t j = 1; j <= reads_len[tasklet_id]; j++) {
						int prior;
						if (tasklet_id == 0 && i == haplotypes_len[haplotype_idx] && j == reads_len[tasklet_id]) {
							result[0] = reads_len[tasklet_id];// (int)reads_arr[tasklet_id][j - 1];
							result[1] = haplotypes_len[haplotype_idx];//(int)haplotypes_arr[haplotype_idx][i - 1];
						}
						if (reads_arr[tasklet_id][j - 1] == haplotypes_arr[haplotype_idx][i - 1]) {
							prior = prior_arr[tasklet_id][2 * (j - 1)];
						}
						else {
							prior = prior_arr[tasklet_id][2 * (j - 1) + 1];
						}
						uint32_t indI = i % MATRIX_LINES;
						uint32_t indI0 = (i - 1) % MATRIX_LINES;

						
						MATCH_CACHE[tasklet_id][indI][j] = fixedAdd( prior , log10SumLog10( log10SumLog10( fixedAdd(MATCH_CACHE[tasklet_id][indI0][j - 1] , transition[matchToMatch]) ,
							fixedAdd(INSERTION_CACHE[tasklet_id][indI0][j - 1] , transition[indelToMatch]) ) , 
							fixedAdd( DELETION_CACHE[tasklet_id][indI0][j - 1], transition[indelToMatch]) ) );

						INSERTION_CACHE[tasklet_id][indI][j] = log10SumLog10( fixedAdd(MATCH_CACHE[tasklet_id][indI][j - 1], transition[matchToInsertion]) , fixedAdd(INSERTION_CACHE[tasklet_id][indI][j - 1], transition[insertionToInsertion]) ) ;

						
						DELETION_CACHE[tasklet_id][indI][j] = log10SumLog10(fixedAdd(MATCH_CACHE[tasklet_id][indI0][j], transition[matchToDeletion]) , fixedAdd(DELETION_CACHE[tasklet_id][indI0][j] , transition[deletionToDeletion]) );
						
					}


					res[tasklet_id] = log10SumLog10( res[tasklet_id], log10SumLog10(MATCH_CACHE[tasklet_id][i % MATRIX_LINES][reads_len[tasklet_id] - 1] , INSERTION_CACHE[tasklet_id][i % MATRIX_LINES][reads_len[tasklet_id] - 1]));
					if (tasklet_id == 0) {
						result[5] = 1;
					}
					//MATCH_CACHE[tasklet_id][0][0] = 0; //needed because we use only two lines of Matrix (and at initialization M[0][0] != M[i][0] i>0)
				}
				if (tasklet_id == 0) {
					for (uint32_t i = 0; i <= haplotypes_len[haplotype_idx]; i++) {
						for (uint32_t j = 0; j<= reads_len[tasklet_id]; j++) {
							mram_write(&MATCH_CACHE[0][i][j], &likelihoods[i][j], sizeof(int));
						}
					}
				}
				//result[2] = result[0] / result[1];
				//result[6] = result[5] / result[1];
				if (tasklet_id == 0) {
					result[6] = 1;
				}

				if (tasklet_id == 0) {
					//mram_write(&res, &likelihoods[haplotype_idx][read_idx], sizeof(res));//TODO: change to write to an array and only task0 writes the whole array to the mram
				}
				if (tasklet_id == 0) {
					result[7] = 1;
				}
				barrier_wait(&my_barrier);
				
			}

		}
		barrier_wait(&my_barrier);
	}
	return 1;
}