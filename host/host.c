#include <dpu.h>
#include <dpu_types.h>
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <unistd.h> 

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

int main() {
	struct dpu_set_t set, dpu;
	uint32_t nr_dpus, each_dpu;
	int result[24];

	FILE* data_file = read_data("1region2.csv");


	sleep(30);

	DPU_ASSERT(dpu_alloc(1, "backend=simulator", &set));
	DPU_ASSERT(dpu_load(set, DPU_BINARY, NULL));
	DPU_ASSERT(dpu_get_nr_dpus(set, &nr_dpus));

	//for (int i = 0; i < nr_regions; i++) {

		populate_mram(set, nr_dpus, 0);

		printf("Launch DPU\n");
		dpu_launch(set, DPU_SYNCHRONOUS);
		printf("Finished DPU work\n");
		DPU_FOREACH(set, dpu, each_dpu) {
			DPU_ASSERT(dpu_copy_from(dpu, "result", 0, &result, sizeof(result)));
			DPU_ASSERT(dpu_copy_from(dpu, "likelihoods", 0, &likelihoods, sizeof(likelihoods)));
			//DPU_ASSERT(dpu_copy_from(dpu, "read_check", 0, &read_check, sizeof(read_check)));
			//DPU_ASSERT(dpu_copy_from(dpu, "read_check2", 0, &read_check2, sizeof(read_check2)));
			//DPU_ASSERT(dpu_copy_from(dpu, "prior_check", 0, &prior_check, sizeof(prior_check)));
			//DPU_ASSERT(dpu_copy_from(dpu, "prior_check2", 0, &prior_check2, sizeof(prior_check2)));
			//printf("[%u] Result at the end: %d and %d %d %d %d %d %d %d\n", each_dpu, result[0], result[1],result[2], result[3], result[4], result[5], result[6], result[7]);
		}
	//}
		
		printf("Likelihhoods\n");
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j <8; j++) {
				printf("%f ", (double)likelihoods[i][j] /ONE);
			}
			printf("\n");
		}
		
		printf("\nResult check:\n");
		for (int i = 0; i < 24; i++) {
			printf("%d ", result[i]);
		}	
		printf("\n");

		/*
		printf("\nREad check2:\n");
		for (int i = 0; i < 120; i++) {
			printf("%f ", read_check2[i]);
		}*/ 
		/*
		printf("\nprior check:\n");
		for (int i = 0; i < 238; i++) {
			printf("%f ", prior_check[i]);
		}
		printf("\nprior check2:\n");
		for (int i = 0; i < 238; i++) {
			printf("%f ", prior_check2[i]);
		}*/

	free_mem(data_file);
	DPU_ASSERT(dpu_free(set));
	return 0; 
}