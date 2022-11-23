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
#include "log.h"
#include "launch_dpus.h"




int64_t* likelihoods[NUMBER_DPUS];

uint64_t nb_cycles[NUMBER_DPUS];

//activate dpus
uint64_t dpu_inactive[NUMBER_DPUS];

//DATA in order to print results
extern uint64_t nr_haplotypes[NR_REGIONS]; //an array keeping number of haplotypes in all regions
extern uint64_t nr_reads[NR_REGIONS]; //idem as haplotypes
extern uint32_t haplotype_region_starts[NUMBER_DPUS][MAX_REGIONS_PER_DPU + 1];
extern uint32_t read_region_starts[NUMBER_DPUS][MAX_REGIONS_PER_DPU + 1];

void init(uint32_t nr_dpus) {
	assert(nr_dpus <= NUMBER_DPUS);
	LOG_DEBUG("allocating likelihood tables (%u tables of size %lu)\n", nr_dpus, sizeof(MAX_HAPLOTYPE_NUM * MAX_READ_NUM));
	for (int i = 0; i < nr_dpus; i++) {
		likelihoods[i] = malloc(MAX_HAPLOTYPE_NUM * MAX_READ_NUM * sizeof(int64_t));
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

	// Start-up parser thread.
	// TODO: implement this

	// Start-up dpu threads.
	LOG_INFO("Launching dpu dispatching\n");
	launch_all_ranks();

	// Start-up result processing thread.
	// TODO: implement this 
	//       if result processing is single-threaded, just call a function without creating a new thread.
	//       So that we don't need to do a sync before returning from main.

	// TODO: Ensure all threads have returned or ensure the result queue is properly closed before this call.
	free_dpus();
	return 0;
}
