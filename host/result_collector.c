#include "buffers.h"
#include "constants.h"

#define MAX_SUB_REGIONS 600
/*
* Collect results from dpu_result_queue
* The results are written in the dpu_result_buffer
* Write results to output file
*/


bool received_regions[TOTAL_REGIONS]; //keeps the regions received
int sub_region_count[TOTAL_REGIONS]; //keeps the number of subregions received - the number of regions
bool allocated_regions[TOTAL_REGIONS]; //bool whether the region was allocated in results
int** results[TOTAL_REGIONS];


void collect_result(FILE* output_file) {
	//initialization
	for (int i = 0; i < TOTAL_REGIONS; i++) {
		received_regions[i] = 0;
		sub_region_count[i] = 0;
		allocated_regions[i] = 0;
	}
	int current_result = queue_take(&dpu_results_queue);
	while (current_result != -1) {
		struct dpu_results_t dpu_results = dpu_results_buffer[current_result];
		for (int i = 0; i < dpu_results.nr_regions; i++) {
			struct region_shape_t region = dpu_results.region_shapes[i];
			//allocate results for region if needed
			if (!allocated_regions[region.region_index]) {
				results[region.region_index] = malloc(region.total_reads_region * sizeof(int*));
				if (!results[region.region_index]) { printf("Error allocating\n"); }
				for (int j = 0; j < region.total_reads_region; j++) {
					results[region.region_index][j] = malloc(region.total_haps_region * sizeof(int));
					if (!results[region.region_index][j]) { printf("Error allocating\n"); }
				}
				allocated_regions[region.region_index] = 1;
			}

			for (int n = 0; n < region.nr_reads; n++) {
				for (int m = 0; m < region.nr_haplotypes; m++) {
					results[region.region_index][n + region.read_offset][m + region.hapl_offset] = dpu_results.likelihoods[n][m];
				}
			}
			sub_region_count[region.region_index]--;
			if (region.last_subregion) {
				sub_region_count[region.region_index] += region.total_nr_subregions;
			}
		}
		queue_release(&dpu_results_queue, current_result);
		current_result = queue_take(&dpu_results_queue);
	}

}

