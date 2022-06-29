#include <dpu.h>
#include <dpu_types.h>
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <unistd.h> 
#include <pthread.h>
#include <semaphore.h>

#include "parser.h"
#include "populateMRAM.h"
#include "constants.h"

#ifndef DPU_BINARY
#define DPU_BINARY "./build/haplotype_dpu"
#endif

#define MAX_HAPLOTYPES 100
#define MAX_READS 250
int likelihoods[MAX_HAPLOTYPES][MAX_READS];

double read_check[50];
double read_check2[120];
double prior_check[238];
double prior_check2[238];

extern uint32_t nr_regions; //number of regions

extern uint32_t* nr_haplotypes; //an array keeping number of haplotypes in all regions
extern uint32_t* nr_reads; //idem as haplotypes

extern uint32_t** reads_len;
extern char*** reads_array;
extern uint32_t*** qualities;
extern uint32_t** haplotypes_len;
extern char*** haplotypes_array;

#define NB_RANKS_MAX 8

typedef struct {
	unsigned int nb_dpus_per_rank[NB_RANKS_MAX];
	unsigned int rank_mram_offset[NB_RANKS_MAX];
	unsigned int nb_dpus;
	struct dpu_set_t all_ranks;
	struct dpu_set_t ranks[NB_RANKS_MAX];
	pthread_mutex_t log_mutex;
	FILE* log_file;
	struct triplet* dpus;
} devices_t;
static devices_t devices;


int main() {
	struct dpu_set_t set, dpu;
	uint32_t nr_dpus, each_dpu;
	int result[24];
	uint32_t nb_cycles, it_counter;
	printf("Will read the data\n");
	FILE* data_file = read_data("1region.csv");
		
	struct dpu_set_t rank;
	unsigned int each_rank;

	DPU_ASSERT(dpu_alloc(DPU_ALLOCATE_ALL, "backend=simulator", &set));

	DPU_RANK_FOREACH(set, rank, each_rank) {
		
		printf("Each rank:%d\n", each_rank);
	}
	DPU_FOREACH(set, dpu, each_dpu) {
		printf("Each dpu=%d\n", each_dpu);
	}

	DPU_ASSERT(dpu_load(set, DPU_BINARY, NULL));
	DPU_ASSERT(dpu_get_nr_dpus(set, &nr_dpus));
	uint32_t a;
	dpu_get_nr_ranks(set, &a);
	printf("Nr_regions is:%d and nr_ranks=%d, nr_dpus=%d \n", nr_regions, a, nr_dpus);
	for (int i = 0; i < nr_regions; i++) {

		populate_mram(set, nr_dpus, i);

		printf("Launch DPU\n");
		dpu_launch(set, DPU_SYNCHRONOUS);
		printf("Finished DPU work\n");

		DPU_FOREACH(set, dpu, each_dpu) {
			DPU_ASSERT(dpu_copy_from(dpu, "result", 0, &result, sizeof(result)));
			DPU_ASSERT(dpu_copy_from(dpu, "likelihoods", 0, &likelihoods, sizeof(likelihoods)));
			DPU_ASSERT(dpu_copy_from(dpu, "nb_cycles", 0, &nb_cycles, sizeof(nb_cycles)));
			DPU_ASSERT(dpu_copy_from(dpu, "it_counter", 0, &it_counter, sizeof(it_counter)));
		}
		printf("Number of cycles: %u for %u reads\n", nb_cycles, it_counter);
		printf("Likelihhoods\n");
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j <18; j++) {
				printf("%f ", (double)likelihoods[i][j] /ONE);
			}
			printf("\n");
		}
		
		printf("\nResult check:\n");
		for (int i = 0; i < 24; i++) {
			printf("%d ", result[i]);
		}	
		printf("\n");
	}

	free_mem(data_file);
	DPU_ASSERT(dpu_free(set));
	return 0; 
}