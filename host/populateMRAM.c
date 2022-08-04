#include "populateMRAM.h"
#include "constants.h"
#include <math.h>



extern uint32_t nr_haplotypes[NR_REGIONS];
extern uint32_t nr_reads[NR_REGIONS];

extern uint32_t offset[NR_REGIONS][OFFSET_SIZE];

extern uint32_t reads_len[TOTAL_READS]; 
extern char reads_array[TOTAL_READS * MAX_READ_LENGTH]; 
extern uint32_t haplotypes_len[TOTAL_HAPS];
extern uint32_t haplotypes_val[TOTAL_HAPS];
extern char haplotypes_array[TOTAL_HAPS * MAX_HAPLOTYPE_LENGTH];
extern uint32_t priors[TOTAL_READS * MAX_READ_LENGTH * 2];



/** 
* populate mram of all the dpus in the given set
* @param iteration is the  
**/
void populate_mram(struct dpu_set_t set, uint32_t nr_dpus, int iteration) {

	struct dpu_set_t dpu; 
	uint32_t each_dpu;
	
	//transfer READS_LEN to DPUs
	DPU_FOREACH(set, dpu, each_dpu) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, &reads_len[ offset[each_dpu][READS_LEN_ARRAY] ]));
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "mram_reads_len", 0, MAX_READ_NUM * sizeof(uint32_t), DPU_XFER_DEFAULT));


	//transfer HAPLOTYPES_LEN to DPUs
	DPU_FOREACH(set, dpu, each_dpu) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, &haplotypes_len[ offset[each_dpu][HAPLOTYPES_LEN_VAL_ARRAY] ]));
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "mram_haplotypes_len", 0, MAX_HAPLOTYPE_NUM * sizeof(uint32_t), DPU_XFER_DEFAULT));

	
	//transfer READS_ARRAY to DPUs
	DPU_FOREACH(set, dpu, each_dpu) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, &reads_array[ offset[each_dpu][READS_ARR] ]));
	}
	uint32_t region_read_size = MAX_READ_NUM * MAX_READ_LENGTH * sizeof(char);
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "mram_reads_array", 0, region_read_size, DPU_XFER_DEFAULT));

	/*
	for (int region = 0; region < NR_REGIONS; region++) {
		for (int i = 0; i < nr_reads[region]; i++) {
			printf("Read New\n");
			printf("%s\n", reads_array[region][i]);
			for (int j = 0; j < reads_len[region][i]; j++) {
				printf("%d | %d | ", priors[region][i][2 * j], priors[region][i][2 * j + 1]);
			}
			printf("\n");
		}
	}*/


	//transfer prior array to each dpu
	//transfer PRIORS to DPUs
	uint32_t prior_read_size = 2 * MAX_READ_NUM * MAX_READ_LENGTH * sizeof(uint32_t);
	DPU_FOREACH(set, dpu, each_dpu) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, &priors[ offset[each_dpu][PRIOR_ARR] ]));
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "mram_priors", 0, prior_read_size, DPU_XFER_DEFAULT));

	//transfer READS_ARRAY to DPUs
	DPU_FOREACH(set, dpu, each_dpu) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, &haplotypes_array[ offset[each_dpu][HAPS_ARR] ]));
	}
	uint32_t haplotypes_arr_size = MAX_HAPLOTYPE_NUM * MAX_HAPLOTYPE_LENGTH * sizeof(char);
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "mram_haplotypes_array",  0, haplotypes_arr_size, DPU_XFER_DEFAULT));


	//transfer HAPLOTYPES_VAL to DPUs
	DPU_FOREACH(set, dpu, each_dpu) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, &haplotypes_val[ offset[each_dpu][HAPLOTYPES_LEN_VAL_ARRAY] ]));
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "mram_haplotypes_val", 0, MAX_HAPLOTYPE_NUM * sizeof(uint32_t), DPU_XFER_DEFAULT));
	

	//transfer NR_HAPLOTYPES to DPUs
	DPU_FOREACH(set, dpu, each_dpu) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, &nr_haplotypes[each_dpu]));
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "nr_haplotypes", 0, sizeof(uint32_t), DPU_XFER_DEFAULT));
	
	//transfer NR_READS to DPUs
	DPU_FOREACH(set, dpu, each_dpu) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, &nr_reads[each_dpu * iteration]));
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "nr_reads", 0, sizeof(uint32_t), DPU_XFER_DEFAULT));



	return;
}
