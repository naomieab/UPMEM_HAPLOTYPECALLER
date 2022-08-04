#include <dpu.h>
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <unistd.h> 
#include <math.h>
#include <time.h>

#include "parser.h"
#include "populateMRAM.h"   
#include "constants.h"

#ifndef DPU_BINARY
#define DPU_BINARY "./build/haplotype_dpu"
#endif

//RESULTS
int likelihoods[NR_REGIONS][MAX_HAPLOTYPE_NUM][MAX_READ_NUM];
int nb_cycles[NR_REGIONS];


//DATA in order to print results
extern uint32_t nr_haplotypes[NR_REGIONS]; //an array keeping number of haplotypes in all regions
extern uint32_t nr_reads[NR_REGIONS]; //idem as haplotypes


 
int main() {
	struct dpu_set_t set, dpu;
	uint32_t nr_dpus, nr_ranks, each_dpu;
	FILE* data_file = read_data("./../HaplotypeCaller/ActiveRegionsReadsHaplotypes.csv");
  time_t start, end;


	DPU_ASSERT(dpu_alloc(DPU_ALLOCATE_ALL, NULL, &set));
	DPU_ASSERT(dpu_load(set, DPU_BINARY, NULL));
	DPU_ASSERT(dpu_get_nr_dpus(set, &nr_dpus));
  DPU_ASSERT(dpu_get_nr_ranks(set, &nr_ranks));
  printf("Nr ranks =%d, and Nr dpus = %d\n", nr_ranks, nr_dpus);
  

	//i is the iteration: if we have several rounds to process on a set of dpus, each iteration process a single round
	for (int iteration = 1; iteration < (NR_REGIONS / nr_dpus) + 1; iteration++) {
    time(&start);
		populate_mram(set, nr_dpus, iteration);
    time(&end);
    printf("Populate MRAM in %ld\n", end-start);
    
		printf("Launch DPU iteration %d\n", iteration);
    time(&start);
		DPU_ASSERT(dpu_launch(set, DPU_SYNCHRONOUS));
    time(&end);
		printf("Finished DPU work: time required %ld\n", end-start);
		 
		DPU_FOREACH(set, dpu, each_dpu) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, &likelihoods[each_dpu * iteration]));
		}
		DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_FROM_DPU, "likelihoods", 0, MAX_HAPLOTYPE_NUM * MAX_READ_NUM * sizeof(uint32_t), DPU_XFER_DEFAULT));


		DPU_FOREACH(set, dpu, each_dpu) {
				DPU_ASSERT(dpu_prepare_xfer(dpu, &nb_cycles[each_dpu * iteration]));
		}
		DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_FROM_DPU, "nb_cycles", 0, sizeof(int), DPU_XFER_DEFAULT));

    /*for (int i = 0; i < nr_dpus; i++) {
			printf("\nREGION %d\n", i*iteration);
			printf("Number of cycles = %d\n", nb_cycles[i*iteration]);
			printf("Number of haplotypes = %d\n", nr_haplotypes[i*iteration]);
			printf("Number of reads = %d\n", nr_reads[i*iteration]);
			for (int j = 0; j < nr_haplotypes[i*iteration]; j++) {
				for (int k = 0; k < nr_reads[i*iteration]; k++) {
					printf("%d => %f | ", likelihoods[i*iteration][j][k], (double)likelihoods[i*iteration][j][k] / (double)ONE);
				}
				printf("\n");
			}
			printf("\n\n\n");
		}*/
		
	}

	free_mem(data_file);
	DPU_ASSERT(dpu_free(set));
	return 0; 
}
