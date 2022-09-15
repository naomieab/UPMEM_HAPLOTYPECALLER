#include "populateMRAM.h"
#include "constants.h"
#include <math.h>



extern uint64_t nr_haplotypes[NUMBER_DPUS];
extern uint64_t nr_reads[NUMBER_DPUS];
extern uint64_t dpu_inactive[NUMBER_DPUS];

extern uint32_t offset[NR_REGIONS][OFFSET_SIZE];
extern uint32_t dpu_region_start_index[NUMBER_DPUS];
extern uint32_t haplotype_region_starts[NUMBER_DPUS][MAX_REGIONS_PER_DPU+1];
extern uint32_t read_region_starts[NUMBER_DPUS][MAX_REGIONS_PER_DPU+1];

extern uint64_t reads_len[TOTAL_READS]; 
extern char reads_array[TOTAL_READS * MAX_READ_LENGTH]; 

extern uint64_t haplotypes_len[TOTAL_HAPS];
extern uint64_t haplotypes_val[TOTAL_HAPS];

extern char haplotypes_array[TOTAL_HAPS * MAX_HAPLOTYPE_LENGTH];
extern uint32_t priors[TOTAL_READS * MAX_READ_LENGTH * 2];
extern int32_t matchToIndel[TOTAL_READS * MAX_READ_LENGTH];


/** 
* populate mram of all the dpus in the given set
* @param iteration is the  
**/
void populate_mram(struct dpu_set_t set, uint32_t nr_dpus, int iteration) {
	struct dpu_set_t dpu; 
	uint32_t each_dpu;
	
	//transfer READS_LEN to DPUs
	DPU_FOREACH(set, dpu, each_dpu) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, &reads_len[ offset[dpu_region_start_index[each_dpu]][READS_LEN_ARRAY] ]));
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "mram_reads_len", 0, MAX_READ_NUM * sizeof(uint64_t), DPU_XFER_DEFAULT));


	//transfer HAPLOTYPES_LEN to DPUs
	DPU_FOREACH(set, dpu, each_dpu) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, &haplotypes_len[ offset[dpu_region_start_index[each_dpu]][HAPLOTYPES_LEN_VAL_ARRAY] ]));
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "mram_haplotypes_len", 0, MAX_HAPLOTYPE_NUM * sizeof(uint64_t), DPU_XFER_DEFAULT));

	
	//transfer READS_ARRAY to DPUs
	DPU_FOREACH(set, dpu, each_dpu) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, &reads_array[ offset[dpu_region_start_index[each_dpu]][READS_ARR] ]));
	}
	uint32_t region_read_size = MAX_READ_NUM * MAX_READ_LENGTH * sizeof(char);
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "mram_reads_array", 0, region_read_size, DPU_XFER_DEFAULT));

	//transfer transitions quals to DPUs
	DPU_FOREACH(set, dpu, each_dpu) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, &matchToIndel[ offset[each_dpu][READS_ARR] ]));
	}
	region_read_size = MAX_READ_NUM * MAX_READ_LENGTH * sizeof(int32_t);
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "mram_matchToIndelArray", 0, region_read_size, DPU_XFER_DEFAULT));


	//transfer prior array to each dpu
	//transfer PRIORS to DPUs
	uint32_t prior_read_size = 2 * MAX_READ_NUM * MAX_READ_LENGTH * sizeof(int32_t);
	DPU_FOREACH(set, dpu, each_dpu) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, &priors[ offset[dpu_region_start_index[each_dpu]][PRIOR_ARR] ]));
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "mram_priors", 0, prior_read_size, DPU_XFER_DEFAULT));

	//transfer HAPLOTYPES_ARRAY to DPUs
	DPU_FOREACH(set, dpu, each_dpu) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, &haplotypes_array[ offset[dpu_region_start_index[each_dpu]][HAPS_ARR] ]));
	}
	uint32_t haplotypes_arr_size = MAX_HAPLOTYPE_NUM * MAX_HAPLOTYPE_LENGTH * sizeof(char);
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "mram_haplotypes_array",  0, haplotypes_arr_size, DPU_XFER_DEFAULT));


	//transfer HAPLOTYPES_VAL to DPUs
	DPU_FOREACH(set, dpu, each_dpu) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, &haplotypes_val[ offset[dpu_region_start_index[each_dpu]][HAPLOTYPES_LEN_VAL_ARRAY] ]));
	}

	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "mram_haplotypes_val", 0, MAX_HAPLOTYPE_NUM * sizeof(uint64_t), DPU_XFER_DEFAULT));

	

	//transfer NR_HAPLOTYPES to DPUs
	DPU_FOREACH(set, dpu, each_dpu) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, &nr_haplotypes[each_dpu]));
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "nr_haplotypes", 0, sizeof(uint64_t), DPU_XFER_DEFAULT));
	
	//transfer NR_READS to DPUs
	DPU_FOREACH(set, dpu, each_dpu) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, &nr_reads[each_dpu]));
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "nr_reads", 0, sizeof(uint64_t), DPU_XFER_DEFAULT));

	//transfer READ_REGION_START_INDEX to DPUs
	DPU_FOREACH(set, dpu, each_dpu) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, &read_region_starts[each_dpu]));
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "read_region_starts", 0, sizeof(uint32_t)*(MAX_REGIONS_PER_DPU+1), DPU_XFER_DEFAULT));

	//transfer HAPLOTYPE_REGION_START_INDEX to DPUs
	DPU_FOREACH(set, dpu, each_dpu) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, &haplotype_region_starts[each_dpu]));
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "haplotype_region_starts", 0, sizeof(uint32_t)*(MAX_REGIONS_PER_DPU+1), DPU_XFER_DEFAULT));

  //transfer dpu_inactive to DPUs
	DPU_FOREACH(set, dpu, each_dpu) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, &dpu_inactive[each_dpu]));
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "dpu_inactive", 0, sizeof(uint64_t), DPU_XFER_DEFAULT));

	return;
}
