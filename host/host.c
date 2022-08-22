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
int64_t likelihoods[NR_REGIONS][MAX_HAPLOTYPE_NUM][MAX_READ_NUM];
uint64_t nb_cycles[NR_REGIONS];


//DATA in order to print results
extern uint64_t nr_haplotypes[NR_REGIONS]; //an array keeping number of haplotypes in all regions
extern uint64_t nr_reads[NR_REGIONS]; //idem as haplotypes



int main() {
	struct dpu_set_t set, dpu;
	uint32_t nr_dpus, nr_ranks, each_dpu;
	FILE* data_file = read_data("./ActiveRegionsReadsHaplotypes.csv");
	time_t start, end;


	DPU_ASSERT(dpu_alloc(DPU_ALLOCATE_ALL, NULL, &set));
	DPU_ASSERT(dpu_load(set, DPU_BINARY, NULL));
	DPU_ASSERT(dpu_get_nr_dpus(set, &nr_dpus));
	DPU_ASSERT(dpu_get_nr_ranks(set, &nr_ranks));
	printf("Nr ranks =%d, and Nr dpus = %d\n", nr_ranks, nr_dpus);


	//i is the iteration: if we have several rounds to process on a set of dpus, each iteration process a single round
	for (int iteration = 0; iteration < (NR_REGIONS / nr_dpus); iteration++) {
		time(&start);
		populate_mram(set, nr_dpus, iteration);
		time(&end);
		printf("Populate MRAM in %ld\n", end - start);

		printf("Launch DPU iteration %d\n", iteration);
		time(&start);
		DPU_ASSERT(dpu_launch(set, DPU_SYNCHRONOUS));
		time(&end);
		printf("Finished DPU work: time required %ld\n", end - start);


		DPU_FOREACH(set, dpu, each_dpu) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, &likelihoods[each_dpu + nr_dpus * iteration]));
		}
		DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_FROM_DPU, "likelihoods", 0, MAX_HAPLOTYPE_NUM * MAX_READ_NUM * sizeof(int64_t), DPU_XFER_DEFAULT));


		DPU_FOREACH(set, dpu, each_dpu) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, &nb_cycles[each_dpu]));
		}
		DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_FROM_DPU, "nb_cycles", 0, sizeof(uint64_t), DPU_XFER_DEFAULT));


		/*for (int i = 0; i < nr_dpus; i++) {
			int region_index = iteration*nr_dpus+i;
			printf("\nREGION %d\n", region_index);
			printf("Number of cycles = %d\n", nb_cycles[region_index]);
			if (region_index < NR_REGIONS) {
				printf("Number of haplotypes = %d\n", nr_haplotypes[region_index]);
				printf("Number of reads = %d\n", nr_reads[region_index]);
				for (int j = 0; j < nr_haplotypes[region_index]; j++) {
					for (int k = 0; k < nr_reads[region_index]; k++) {
						printf("%d => %f | ", likelihoods[region_index][j][k], (double)likelihoods[region_index][j][k] / (double)ONE);
					}
					printf("\n");
				}
				printf("\n\n\n");
			}
		}*/

	}

	free_mem(data_file);
	DPU_ASSERT(dpu_free(set));
	return 0;
}
