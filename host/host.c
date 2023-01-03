#include <dpu.h>
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <unistd.h> 
#include <math.h>
#include <time.h>

#include "buffers.h"
#include "parser.h"
#include "constants.h"
#include "log.h"
#include "launch_dpus.h"
#include "result_collector.h"



int64_t* likelihoods[NUMBER_DPUS];

uint64_t nb_cycles[NUMBER_DPUS];

//activate dpus
uint64_t dpu_inactive[NUMBER_DPUS];

//DATA in order to print results
extern uint64_t nr_haplotypes[NR_REGIONS]; //an array keeping number of haplotypes in all regions
extern uint64_t nr_reads[NR_REGIONS]; //idem as haplotypes
extern uint32_t haplotype_region_starts[NUMBER_DPUS][MAX_REGIONS_PER_DPU + 1];
extern uint32_t read_region_starts[NUMBER_DPUS][MAX_REGIONS_PER_DPU + 1];

// TODO: This temporary function has been used for testing. Delete it once it is not useful anymore
void* send_dummy_region() {
	reads_len_buffer[0] = 5;
	char read[10] = "ACTGC";
	char hapl[10] = "ACAGC";
	memcpy(reads_array_buffer, read, 6);
	haplotypes_len_buffer[0] = 5;
	haplotypes_val_buffer[0] = 1;
	memcpy(haplotypes_array_buffer, hapl, 6);
	int index = queue_put(&dpu_regions_queue);
	dpu_regions_buffer[index].dpu_inactive = false;
	dpu_regions_buffer[index].first_region_index = 0;
	dpu_regions_buffer[index].nr_regions = 1;
	dpu_regions_buffer[index].nr_reads = 1;
	dpu_regions_buffer[index].nr_haplotypes = 1;
	dpu_regions_buffer[index].region_shapes[0].region_index = 0;
	dpu_regions_buffer[index].region_shapes[0].nr_reads = 1;
	dpu_regions_buffer[index].region_shapes[0].nr_haplotypes = 1;
	dpu_regions_buffer[index].region_shapes[0].read_offset = 0;
	dpu_regions_buffer[index].region_shapes[0].hapl_offset = 0;
	dpu_regions_buffer[index].region_shapes[0].total_nr_subregions = 1;
	dpu_regions_buffer[index].haplotype_region_starts[0] = 0;
	dpu_regions_buffer[index].haplotype_region_starts[1] = 1;
	dpu_regions_buffer[index].read_region_starts[0] = 0;
	dpu_regions_buffer[index].read_region_starts[1] = 1;
	dpu_regions_buffer[index].reads_len = reads_len_buffer;
	dpu_regions_buffer[index].reads_array = reads_array_buffer;
	dpu_regions_buffer[index].haplotypes_len = haplotypes_len_buffer;
	dpu_regions_buffer[index].haplotypes_val = haplotypes_val_buffer;
	dpu_regions_buffer[index].haplotypes_array = haplotypes_array_buffer;
	dpu_regions_buffer[index].priors = priors;
	dpu_regions_buffer[index].match_to_indel = match_to_indel_buffer;

	queue_make_available(&dpu_regions_queue, index);
	queue_close(&dpu_regions_queue, MAX_RANKS);
 return NULL;
}

// TODO: This temporary function has been used for testing. Delete it once it is not useful anymore
void print_all_dpu_results() {
	while (true) {
		int index = queue_take(&dpu_results_queue);
		if (index < 0) {
			LOG_INFO("results queue closed\n");
			return;
		}
		printf("got result from dpus: %d regions computed in %u cycles\n", dpu_results_buffer[index].nr_regions, dpu_results_buffer[index].nr_cycles);
		queue_release(&dpu_results_queue, index);
	}
}

int main(int argc, char* argv[]) {
	if (argc != 4) {
		LOG_FATAL("Wrong parameters! You must provide:\n-1st argument: INPUT regions file to read (sorted)\n-2nd argument: OUTPUT file to write results computed by DPUs\n-3rd argument: PERF file to write DPUs performances (csv file)\nThe line should be make test INPUT=input_file.csv OUTPUT=output_file.txt PERF=performance_file.csv\n");
		return 0;
	}
	FILE* csv_result = fopen(argv[3], "w");
	if (!csv_result) {
		LOG_FATAL("Can't open result file\n");
		return 0;
	}
	FILE* data_file = fopen(argv[1], "r");
	if (!data_file) {
		LOG_FATAL("Can't read input file");
		return 0;
	}
	FILE* result_file = fopen(argv[2], "w");
	if (!result_file) {
		LOG_FATAL("Can't open result file");
		return 0;
	}
	//Set up queues:
	queue_init(&dpu_regions_queue, DPU_INPUT_BUFFER_SIZE);
	queue_init(&dpu_results_queue, DPU_OUTPUT_BUFFER_SIZE);

	// Start-up parser thread.
	// TODO: implement this
	LOG_INFO("Launching parser\n");
	pthread_t parser_thread;
	pthread_create(&parser_thread, NULL, read_data, (void*)data_file);
	//send_dummy_region();
 //pthread_create(&parser_thread, NULL, send_dummy_region, NULL);

	// Start-up dpu threads.
	LOG_INFO("Launching dpu dispatching\n");
	launch_all_ranks();

	// Start-up result processing thread.
	// TODO: implement this 
	//       if result processing is single-threaded, just call a function without creating a new thread.
	//       So that we don't need to do a sync before returning from main.
	LOG_INFO("Launching result printer\n");
	//print_all_dpu_results();
	collect_result(result_file);

	// TODO: Ensure all threads have returned or ensure the result queue is properly closed before this call.
	free_dpus();
	return 0;
}


