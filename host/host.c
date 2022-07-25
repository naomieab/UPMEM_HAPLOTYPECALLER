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


int likelihoods[NR_REGIONS][MAX_HAPLOTYPE_NUM][MAX_READ_NUM];

int nb_cycles[NR_REGIONS];

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
	FILE* data_file = read_data("ActiveRegionsDetails.csv");


	

	DPU_ASSERT(dpu_alloc(DPU_ALLOCATE_ALL, "backend=simulator", &set));
	DPU_ASSERT(dpu_load(set, DPU_BINARY, NULL));
	DPU_ASSERT(dpu_get_nr_dpus(set, &nr_dpus));
	uint32_t a;
	dpu_get_nr_ranks(set, &a);
	
	//i is the iteration: if we have several rounds to process on a set of dpus, each iteration process a single round
	for (int iteration = 1; iteration < (NR_REGIONS / nr_dpus) + 1; iteration++) {
		populate_mram(set, nr_dpus, iteration);

		printf("Launch DPU\n");
		DPU_ASSERT(dpu_launch(set, DPU_SYNCHRONOUS));
		printf("Finished DPU work\n");

		DPU_FOREACH(set, dpu, each_dpu) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, &likelihoods[each_dpu * iteration]));
		}
		DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_FROM_DPU, "likelihoods", 0, MAX_HAPLOTYPE_NUM * MAX_READ_NUM * sizeof(uint32_t), DPU_XFER_DEFAULT));


		DPU_FOREACH(set, dpu, each_dpu) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, &nb_cycles[each_dpu * iteration]));
		}
		DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_FROM_DPU, "nb_cycles", 0, sizeof(int), DPU_XFER_DEFAULT));

		for (int i = 0; i < NR_REGIONS; i++) {
			printf("Number of cycles = %d\n", nb_cycles[i]);
			printf("Number of reads = %d\n", nr_reads[i]);
			printf("Number of haplotypes = %d\n", nr_haplotypes[i]);
			for (int j = 0; j < nr_haplotypes[i]; j++) {
				for (int k = 0; k < nr_reads[i]; k++) {
					printf("%d | ", likelihoods[i][j][k]);
				}
				printf("\n");
			}
			printf("\n\n\n");
		}
		
		
		
	}

	free_mem(data_file);
	DPU_ASSERT(dpu_free(set));
	return 0; 
}