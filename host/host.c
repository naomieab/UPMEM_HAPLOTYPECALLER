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

//activate dpus
uint64_t dpu_inactive[NR_REGIONS];

//DATA in order to print results
extern uint64_t nr_haplotypes[NR_REGIONS]; //an array keeping number of haplotypes in all regions
extern uint64_t nr_reads[NR_REGIONS]; //idem as haplotypes



int main(int argc, char* argv[]) {
	if (argc != 4) {
		printf("Wrong parameters! You must provide:\n-1st argument: INPUT regions file to read (sorted)\n-2nd argument: OUTPUT file to write results computed by DPUs\n-3rd argument: PERF file to write DPUs performances (csv file)\nThe line should be make test INPUT=input_file.csv OUTPUT=output_file.txt PERF=performance_file.csv\n");
		return 0;
	}
	struct dpu_set_t set, dpu;
	uint32_t nr_dpus, nr_ranks, each_dpu;
	time_t start, end, dpu_time = 0;

	FILE* csv_result = fopen(argv[3], "w");
	if (!csv_result) {
		printf("Can't open result file\n");
		return 0;
	}
	FILE* data_file = fopen(argv[1], "r");
	if (!data_file) {
		printf("Can't read input file");
		return 0;
	}
	FILE* result_file = fopen(argv[2], "w");
	if (!result_file) {
		printf("Can't open result file");
		return 0;
	}


	DPU_ASSERT(dpu_alloc(DPU_ALLOCATE_ALL, NULL, &set));
	DPU_ASSERT(dpu_load(set, DPU_BINARY, NULL));
	DPU_ASSERT(dpu_get_nr_dpus(set, &nr_dpus));
	DPU_ASSERT(dpu_get_nr_ranks(set, &nr_ranks));
	printf("Nr ranks =%d, and Nr dpus = %d\n", nr_ranks, nr_dpus);
	printf("Number of iterations: %d\n", (int)(ceil((double)TOTAL_REGIONS / nr_dpus)));
	int dpu_iterations = (int)(ceil((double)TOTAL_REGIONS / nr_dpus));
	for (int i = 0; i < NR_REGIONS; i++) { dpu_inactive[i] = 0; }

	//i is the iteration: if we have several rounds to process on a set of dpus, each iteration process a single round
	for (int iteration = 0; iteration < dpu_iterations; iteration++) {
		printf("Starting iteration: %d\n", iteration);
		if (iteration == dpu_iterations-1) {
			nr_dpus = TOTAL_REGIONS % NR_REGIONS;
			for (int i = nr_dpus; i < NR_REGIONS; i++) { dpu_inactive[i] = 1; }
		}
		data_file = read_data(data_file, nr_dpus);

		populate_mram(set, nr_dpus, iteration);

		fprintf(stderr, "Launch DPU iteration %d\n", iteration);
		time(&start);
		DPU_ASSERT(dpu_launch(set, DPU_SYNCHRONOUS));
		time(&end);
		fprintf(stderr, "Finished DPU work: time required %ld\n", end - start);
		dpu_time += (end - start);


		DPU_FOREACH(set, dpu, each_dpu) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, &likelihoods[each_dpu]));
		}
		DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_FROM_DPU, "likelihoods", 0, MAX_HAPLOTYPE_NUM * MAX_READ_NUM * sizeof(int64_t), DPU_XFER_DEFAULT));

		DPU_FOREACH(set, dpu, each_dpu) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, &nb_cycles[each_dpu]));
		}
		DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_FROM_DPU, "nb_cycles", 0, sizeof(uint64_t), DPU_XFER_DEFAULT));


		for (int i = 0; i < nr_dpus; i++) {
			int region_index = iteration*nr_dpus+i;
			//fprintf(result_file, "\nREGION %d\n", region_index);
			//fprintf(csv_result, "%d\n", nb_cycles[i]);
			//fprintf(result_file, "Number of haplotypes = %d\n", nr_haplotypes[i]);
			//fprintf(result_file, "Number of reads = %d\n", nr_reads[i]);
			//fprintf(result_file, "First likelihood : %d, %d, %d, %d\n", likelihoods[i][0][0], likelihoods[i][0][1], likelihoods[i][1][0], likelihoods[i][1][1]);
			for (int k = 0; k < nr_reads[i]; k++) {
				for (int j = 0; j < nr_haplotypes[i]; j++) {
					//printf("%d => %f | ", likelihoods[i][j][k], (double)likelihoods[i][j][k] / (double)ONE);
		  			fprintf(result_file, "%f, ", (double)likelihoods[i][j][k] / (double)ONE);
				}
				fprintf(result_file,"\n");
			}
			fprintf(result_file, "\n\n\n");
		}

	}
	fclose(result_file);
	fclose(csv_result);
	free_mem(data_file);
	DPU_ASSERT(dpu_free(set));
	printf("TOTAL DPU TIME: %ld seconds\n", dpu_time);
	return 0;
}
