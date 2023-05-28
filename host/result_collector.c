#include "buffers.h"
#include "constants.h"
#include <stdio.h>
#include <stdlib.h>
#include "log.h"

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
     //LOG_INFO("Writing results %d; %d regions\n", current_result, dpu_results_buffer[current_result].nr_regions);
		int likelihood_hap_idx = 0, likelihood_read_idx = 0;
    for (int i = 0; i < dpu_results_buffer[current_result].nr_regions; i++) {
      //LOG_INFO("Total:%d h %d r, Curr %d h %d r, Offset %d h %d r\n", dpu_results_buffer[current_result].region_shapes[i].total_haps_region, dpu_results_buffer[current_result].region_shapes[i].total_reads_region, dpu_results_buffer[current_result].region_shapes[i].nr_haplotypes,dpu_results_buffer[current_result].region_shapes[i].nr_reads, dpu_results_buffer[current_result].region_shapes[i].hapl_offset, dpu_results_buffer[current_result].region_shapes[i].read_offset);
			//struct region_shape_t region = dpu_results.region_shapes[i];
			//allocate results for region if needed
			// printf("Results of region %d with %d reads and %d haps for region %d, last subregion %d and nb subregions %d\n",dpu_results_buffer[current_result].region_shapes[i].region_index,  dpu_results_buffer[current_result].region_shapes[i].nr_reads, dpu_results_buffer[current_result].region_shapes[i].nr_haplotypes, dpu_results_buffer[current_result].region_shapes[i].region_index,  dpu_results_buffer[current_result].region_shapes[i].last_subregion,  dpu_results_buffer[current_result].region_shapes[i].total_nr_subregions);
      if (!allocated_regions[dpu_results_buffer[current_result].region_shapes[i].region_index]) {
				region_reads_nb[dpu_results_buffer[current_result].region_shapes[i].region_index] = dpu_results_buffer[current_result].region_shapes[i].total_reads_region;
				region_haps_nb[dpu_results_buffer[current_result].region_shapes[i].region_index] = dpu_results_buffer[current_result].region_shapes[i].total_haps_region;
				results[dpu_results_buffer[current_result].region_shapes[i].region_index] = malloc(dpu_results_buffer[current_result].region_shapes[i].total_reads_region * sizeof(int*));
				if (!results[dpu_results_buffer[current_result].region_shapes[i].region_index]) { LOG_ERROR("Error allocating\n"); }
				for (int j = 0; j < dpu_results_buffer[current_result].region_shapes[i].total_reads_region; j++) {
					results[dpu_results_buffer[current_result].region_shapes[i].region_index][j] = malloc(dpu_results_buffer[current_result].region_shapes[i].total_haps_region * sizeof(int));
					if (!results[dpu_results_buffer[current_result].region_shapes[i].region_index][j]) { LOG_ERROR("Error allocating\n"); }
				}
				allocated_regions[dpu_results_buffer[current_result].region_shapes[i].region_index] = 1;
        //LOG_INFO("ALLOCATED \n");
			}

			for (int n = 0; n < dpu_results_buffer[current_result].region_shapes[i].nr_reads; n++) {
				for (int m = 0; m < dpu_results_buffer[current_result].region_shapes[i].nr_haplotypes; m++) {
					results[dpu_results_buffer[current_result].region_shapes[i].region_index][n + dpu_results_buffer[current_result].region_shapes[i].read_offset][m + dpu_results_buffer[current_result].region_shapes[i].hapl_offset] = dpu_results_buffer[current_result].likelihoods[(likelihood_hap_idx+m)*MAX_READ_NUM+(likelihood_read_idx+n)];
				}
			}
      likelihood_hap_idx+=dpu_results_buffer[current_result].region_shapes[i].nr_haplotypes;
      likelihood_read_idx+=dpu_results_buffer[current_result].region_shapes[i].nr_reads;
			sub_region_count[dpu_results_buffer[current_result].region_shapes[i].region_index]--;
			if (dpu_results_buffer[current_result].region_shapes[i].last_subregion) {
				sub_region_count[dpu_results_buffer[current_result].region_shapes[i].region_index] += dpu_results_buffer[current_result].region_shapes[i].total_nr_subregions;
			}
			if (sub_region_count[dpu_results_buffer[current_result].region_shapes[i].region_index] == 0) {
				received_regions[dpu_results_buffer[current_result].region_shapes[i].region_index] = 1;
			}
		}
   
		queue_release(&dpu_results_queue, current_result);

        LOG_DEBUG("queue take results\n");
		current_result = queue_take(&dpu_results_queue);
		while (received_regions[next_to_write] == 1 && next_to_write < TOTAL_REGIONS) {
			//write results to output and free mem
			for (int i = 0; i < region_reads_nb[next_to_write]; i++) {
				for (int j = 0; j < region_haps_nb[next_to_write]; j++) {
					fprintf(output_file, "%f,", (double)results[next_to_write][i][j] / (double)ONE);
				}
				fprintf(output_file, "\n");
				// LOG_DEBUG("Free something\n");
				free(results[next_to_write][i]);
			}
			fprintf(output_file, "\n\n\n");
			free(results[next_to_write]);
			next_to_write++;
		}
	}
}

