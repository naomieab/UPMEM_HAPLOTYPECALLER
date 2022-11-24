#include "buffers.h"
#include "constants.h"

#define MAX_SUB_REGIONS 600
/*
* Collect results from dpu_result_queue
* The results are written in the dpu_result_buffer
* Write results to output file
*/


bool received_regions[TOTAL_REGIONS]; //keeps the regions received
int sub_region_count[TOTAL_REGIONS]; //keeps the number of subregions received
bool allocated_regions[TOTAL_REGIONS]; //bool whether the region was allocated in results
int** results[TOTAL_REGIONS];

void collect_result(FILE* output_file) {
	int current_region = queue_take(&dpu_results_queue);
	while (current_region != -1) {


		current_region = queue_take(&dpu_results_queue);
	}

}

