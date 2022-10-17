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


int64_t *likelihoods[NUMBER_DPUS];
// int64_t likelihoods[NUMBER_DPUS][MAX_HAPLOTYPE_NUM][MAX_READ_NUM];

uint64_t nb_cycles[NUMBER_DPUS];

//activate dpus
uint64_t dpu_inactive[NUMBER_DPUS];

//DATA in order to print results
extern uint64_t nr_haplotypes[NR_REGIONS]; //an array keeping number of haplotypes in all regions
extern uint64_t nr_reads[NR_REGIONS]; //idem as haplotypes
extern uint32_t haplotype_region_starts[NUMBER_DPUS][MAX_REGIONS_PER_DPU+1];
extern uint32_t read_region_starts[NUMBER_DPUS][MAX_REGIONS_PER_DPU+1];

void init(uint32_t nr_dpus) {
	assert(nr_dpus<=NUMBER_DPUS);
	printf("allocating likelihood tables (%u tables of size %lu)\n", nr_dpus, sizeof(MAX_HAPLOTYPE_NUM*MAX_READ_NUM));
	for (int i=0; i<nr_dpus; i++) {
		likelihoods[i] = malloc(MAX_HAPLOTYPE_NUM*MAX_READ_NUM*sizeof(int64_t));
	}
}

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
	int global_region_index = 0;

    unsigned long int total_cycles = 0;

	for (int i = 0; i < NUMBER_DPUS; i++) { dpu_inactive[i] = 0; }

	init(nr_dpus);
	
	//i is the iteration: if we have several rounds to process on a set of dpus, each iteration process a single round
	for (int iteration = 0; iteration < dpu_iterations; iteration++) {
		printf("Starting iteration: %d\n", iteration);
		if (iteration == dpu_iterations-1) {
			nr_dpus = TOTAL_REGIONS % NR_REGIONS;
			for (int i = nr_dpus; i < NUMBER_DPUS; i++) { dpu_inactive[i] = 1; }
		}
		data_file = read_data(data_file, nr_dpus);
        if (nr_reads[0]==0) {
            // If first dpu is empty, then all dpus are which means we have reached the end of the file.
            printf("end of file\n");
            break;
        }

		populate_mram(set, nr_dpus, iteration);

		fprintf(stderr, "Launch DPU iteration %d\n", iteration);
		time(&start);
		DPU_ASSERT(dpu_launch(set, DPU_SYNCHRONOUS));
		time(&end);
		fprintf(stderr, "Finished DPU work: time required %ld\n", end - start);
		dpu_time += (end - start);
		 
		
		DPU_FOREACH(set, dpu, each_dpu) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, likelihoods[each_dpu]));
		}
		DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_FROM_DPU, "likelihoods", 0, MAX_HAPLOTYPE_NUM* MAX_READ_NUM * sizeof(int64_t), DPU_XFER_DEFAULT));

		DPU_FOREACH(set, dpu, each_dpu) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, &nb_cycles[each_dpu]));
		}
		DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_FROM_DPU, "nb_cycles", 0, sizeof(uint64_t), DPU_XFER_DEFAULT));

		fprintf(stderr, "Finished transfering data back\n");
		
		for (int i = 0; i < nr_dpus; i++) {
			int local_region_index = -1;
			// fprintf(result_file, "\nDPU %d\n", i);
			fprintf(csv_result, "%d\n", nb_cycles[i]);
            total_cycles += nb_cycles[i];
			// fprintf(result_file, "Number of haplotypes = %d\n", nr_haplotypes[i]);
			// fprintf(result_file, "Number of reads = %d\n", nr_reads[i]);
			// fprintf(result_file, "First likelihood : %d, %d, %d, %d\n", likelihoods[i][0][0], likelihoods[i][0][1], likelihoods[i][1][0], likelihoods[i][1][1]);
			for (int k = 0; k < nr_reads[i]; k++) {

				if (k >= read_region_starts[i][local_region_index+1]) {
					local_region_index++;
					global_region_index++;
					// fprintf(result_file, "\nREGION %d\n", global_region_index);
					// fprintf(result_file, "Number of haplotypes = %d\n", haplotype_region_starts[i][local_region_index+1]-haplotype_region_starts[i][local_region_index]);
					// fprintf(result_file, "Number of reads = %d\n", read_region_starts[i][local_region_index+1]-read_region_starts[i][local_region_index]);
				}
				for (int j = haplotype_region_starts[i][local_region_index]; j < haplotype_region_starts[i][local_region_index+1]; j++) {
					//printf("%d => %f | ", likelihoods[i][j][k], (double)likelihoods[i][j][k] / (double)ONE);
		  			fprintf(result_file, "%f, ", (double)likelihoods[i][j*MAX_READ_NUM + k] / (double)ONE);

				}
				fprintf(result_file,"\n");
			}
		//	fprintf(result_file, "\n\n\n");
		}

	}
    printf("total cycles: %lu\n", total_cycles);
    fflush(stdout);
	fclose(result_file);
	fclose(csv_result);
	free_mem(data_file);
	DPU_ASSERT(dpu_free(set));
	printf("TOTAL DPU TIME: %ld seconds\n", dpu_time);
	return 0;
}
