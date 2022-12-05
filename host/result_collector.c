#include "buffers.h"
#include "constants.h"
#include <stdio.h>
#include <stdlib.h>


/*
* Collect results from dpu_result_queue
* The results are written in the dpu_result_buffer
* Write results to output file
*/


bool received_regions[TOTAL_REGIONS]; //keeps the regions received
int sub_region_count[TOTAL_REGIONS]; //keeps the number of subregions received - the number of regions
bool allocated_regions[TOTAL_REGIONS]; //bool whether the region was allocated in results
int** results[TOTAL_REGIONS];
int region_reads_nb[TOTAL_REGIONS];
int region_haps_nb[TOTAL_REGIONS];


void collect_result(FILE* output_file) {
	//initialization
	int next_to_write = 0;
	for (int i = 0; i < TOTAL_REGIONS; i++) {
		sub_region_count[i] = 0;
		allocated_regions[i] = 0;
	}
	int current_result = queue_take(&dpu_results_queue);
	while (current_result != -1) {
		for (int i = 0; i < dpu_results_buffer[current_result].nr_regions; i++) {
      
			//struct region_shape_t region = dpu_results.region_shapes[i];
			//allocate results for region if needed
			printf("Results with %d reads and %d haps\n", dpu_results_buffer[current_result].region_shapes[i].total_reads_region, dpu_results_buffer[current_result].region_shapes[i].total_haps_region);
      if (!allocated_regions[dpu_results_buffer[current_result].region_shapes[i].region_index]) {
				region_reads_nb[dpu_results_buffer[current_result].region_shapes[i].region_index] = dpu_results_buffer[current_result].region_shapes[i].total_reads_region;
				region_haps_nb[dpu_results_buffer[current_result].region_shapes[i].region_index] = dpu_results_buffer[current_result].region_shapes[i].total_haps_region;
				results[dpu_results_buffer[current_result].region_shapes[i].region_index] = malloc(dpu_results_buffer[current_result].region_shapes[i].total_reads_region * sizeof(int*));
				if (!results[dpu_results_buffer[current_result].region_shapes[i].region_index]) { printf("Error allocating\n"); }
				for (int j = 0; j < dpu_results_buffer[current_result].region_shapes[i].total_reads_region; j++) {
					results[dpu_results_buffer[current_result].region_shapes[i].region_index][j] = malloc(dpu_results_buffer[current_result].region_shapes[i].total_haps_region * sizeof(int));
					if (!results[dpu_results_buffer[current_result].region_shapes[i].region_index][j]) { printf("Error allocating\n"); }
				}
				allocated_regions[dpu_results_buffer[current_result].region_shapes[i].region_index] = 1;
			}

			for (int n = 0; n < dpu_results_buffer[current_result].region_shapes[i].nr_reads; n++) {
				for (int m = 0; m < dpu_results_buffer[current_result].region_shapes[i].nr_haplotypes; m++) {
					results[dpu_results_buffer[current_result].region_shapes[i].region_index][n + dpu_results_buffer[current_result].region_shapes[i].read_offset][m + dpu_results_buffer[current_result].region_shapes[i].hapl_offset] = dpu_results_buffer[current_result].likelihoods[n][m];
				}
			}
			sub_region_count[dpu_results_buffer[current_result].region_shapes[i].region_index]--;
			if (dpu_results_buffer[current_result].region_shapes[i].last_subregion) {
				sub_region_count[dpu_results_buffer[current_result].region_shapes[i].region_index] += dpu_results_buffer[current_result].region_shapes[i].total_nr_subregions;
			}
			if (sub_region_count[dpu_results_buffer[current_result].region_shapes[i].region_index] == 0) {
				received_regions[dpu_results_buffer[current_result].region_shapes[i].region_index] = 1;
			}
		}
		queue_release(&dpu_results_queue, current_result);
		current_result = queue_take(&dpu_results_queue);
		while (received_regions[next_to_write] == 1 && next_to_write < TOTAL_REGIONS) {
			//write results to output and free mem
			for (int i = 0; i < region_reads_nb[next_to_write]; i++) {
				for (int j = 0; j < region_haps_nb[next_to_write]; j++) {
					fprintf(output_file, "%f,", (double)results[next_to_write][i][j] / (double)ONE);
				}
				fprintf(output_file, "\n");
				printf("Free something\n");
				free(results[next_to_write][i]);
			}
			fprintf(output_file, "\n\n\n");
			free(results[next_to_write]);
			next_to_write++;
		}
	}
}

