#include "populateMRAM.h"
#include "constants.h"

extern uint32_t nr_regions; //number of regions


extern char** reads;
extern char** haplotypes;
extern int* total_read_length;
extern int* total_hap_length;

extern uint32_t nr_haplotypes[NR_REGIONS];
extern uint32_t nr_reads[NR_REGIONS];

extern uint32_t reads_len[NR_REGIONS][MAX_READ_NUM];
extern char reads_array[NR_REGIONS][MAX_READ_NUM][MAX_READ_LENGTH];
extern uint32_t qualities[NR_REGIONS][MAX_READ_NUM][MAX_READ_LENGTH];
extern uint32_t haplotypes_len[NR_REGIONS][MAX_HAPLOTYPE_NUM];
extern uint32_t haplotypes_val[NR_REGIONS][MAX_HAPLOTYPE_NUM];
extern char haplotypes_array[NR_REGIONS][MAX_HAPLOTYPE_NUM][MAX_HAPLOTYPE_LENGTH];

uint32_t priors[NR_REGIONS][MAX_READ_NUM][2 * MAX_READ_LENGTH];

void free_prior(int** array, int region) {
	for (int i = 0; i < nr_reads[region]; i++) {
		free(array[i]);
	}
	free(array);
	return;
}


void print(int** prior, int region) {
	for (int i = 0; i < nr_reads[region]; i++) {
		for (int j = 0; j < reads_len[region][i] - 1; j++) {
			printf("%d %d ", prior[i][2 * j], prior[i][2 * j + 1]);
		}
		printf("\n");
	}
}

/** 
* populate mram of all the dpus in the given set
* @param iteration is the 
**/
void populate_mram(struct dpu_set_t set, uint32_t nr_dpus, int iteration) {

	struct dpu_set_t dpu;
	uint32_t each_dpu;
	uint32_t offset = 0;
	uint32_t offsets[6]; //6 is the number of arrays we need to copy to the mram
	offsets[0] = 0;
	


	DPU_FOREACH(set, dpu, each_dpu) {

		DPU_ASSERT(dpu_prepare_xfer(dpu, &reads_len[each_dpu * iteration]));
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "mram_reads_len", 0, MAX_READ_NUM * sizeof(uint32_t), DPU_XFER_DEFAULT));

	offset += MAX_READ_NUM * sizeof(uint32_t);
	offsets[1] = offset;

	DPU_FOREACH(set, dpu, each_dpu) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, &haplotypes_len[each_dpu * iteration]));
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "mram_haplotypes_len", 0, MAX_HAPLOTYPE_NUM * sizeof(uint32_t), DPU_XFER_DEFAULT));

	offset += MAX_HAPLOTYPE_NUM * sizeof(uint32_t);
	offsets[2] = offset;
	

	int reads_arr_size = 0;
	
	

	DPU_FOREACH(set, dpu, each_dpu) {

		DPU_ASSERT(dpu_prepare_xfer(dpu, &reads_array[each_dpu * iteration]));
	}
	uint32_t region_read_size = MAX_READ_NUM * MAX_READ_LENGTH * sizeof(char);
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "mram_reads_array", 0, region_read_size, DPU_XFER_DEFAULT));
	offset += region_read_size;
	offsets[3] = offset;

	
	
	//create prior array to transfer
	for (int region = 0; region < NR_REGIONS; region++) {
		for (int i = 0; i < MAX_READ_NUM; i++) {
			for (int j = 0; j < reads_len[region][i]; j += 1) {
				double prior = pow((double)10, -(double)qualities[region][i][j] / 10.0);
				priors[region][i][2 * j] = (int)((1 - prior) * ONE);
				priors[region][i][2 * j + 1] = (int)((prior / 3) * ONE);

			}
		}
	}

	//transfer prior array to each dpu
	uint32_t prior_read_size = 2 * MAX_READ_NUM * MAX_READ_LENGTH * sizeof(uint32_t);
	DPU_FOREACH(set, dpu, each_dpu) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, &priors[each_dpu * iteration]));
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "mram_priors", 0, prior_read_size, DPU_XFER_DEFAULT));
	offset += prior_read_size;
	offsets[4] = offset;


	int hap_arr_size = 0;

	DPU_FOREACH(set, dpu, each_dpu) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, &haplotypes_len[each_dpu * iteration]));
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "mram_haplotypes_len", 0, MAX_HAPLOTYPE_NUM * sizeof(uint32_t), DPU_XFER_DEFAULT));
	offset += MAX_HAPLOTYPE_NUM * sizeof(uint32_t);
	offsets[5] = offset;

	DPU_FOREACH(set, dpu, each_dpu) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, &haplotypes_val[each_dpu * iteration]));
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "mram_haplotypes_val", 0, MAX_HAPLOTYPE_NUM * sizeof(uint32_t), DPU_XFER_DEFAULT));
	offset += MAX_HAPLOTYPE_NUM * sizeof(uint32_t);
	offsets[5] = offset;
	

//TODO: ADD NR READS and NR HAPLOTYPES TRANSFER TO DPUS
	DPU_FOREACH(set, dpu, each_dpu) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, &nr_haplotypes[each_dpu * iteration]));
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "nr_haplotypes", 0, sizeof(uint32_t), DPU_XFER_DEFAULT));
	
	DPU_FOREACH(set, dpu, each_dpu) {
		DPU_ASSERT(dpu_prepare_xfer(dpu, &nr_reads[each_dpu * iteration]));
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "nr_reads", 0, sizeof(uint32_t), DPU_XFER_DEFAULT));

	return;
}
