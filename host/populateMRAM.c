#include "populateMRAM.h"
#include "constants.h"

extern uint32_t nr_regions; //number of regions

extern uint32_t* nr_haplotypes; //an array keeping number of haplotypes in all regions
extern uint32_t* nr_reads; //idem as haplotypes

extern uint32_t** reads_len;
extern char*** reads_array;
extern char** reads;
extern uint32_t*** qualities;
extern uint32_t** haplotypes_len;
extern uint32_t** haplotypes_val; 
extern char*** haplotypes_array;
extern char** haplotypes;
extern int* total_read_length;
extern int* total_hap_length;



void free_prior(int** array, int region) {
	for (int i = 0; i < nr_reads[region]; i++) {
		free(array[i]);
	}
	free(array);
	return;
}

int** create_prior(int region) {
	int** prior_array = malloc(nr_reads[region] * sizeof(int*));
	//uint32_t** prior_integer = malloc(nr_reads[region] * sizeof(uint32_t*));
	int read_len;
	for (int i = 0; i < nr_reads[region]; i++) {
		read_len = reads_len[region][i];
		prior_array[i] = malloc(2 * read_len * sizeof(int));
		//prior_integer[i] = malloc(2 * (read_len - 1) * sizeof(uint32_t));
		for (int j = 0; j < read_len; j += 1) {
			double prior = pow((double)10, -(double)qualities[region][i][j] / 10.0);
			//float prior2 = pow((double)10, -(double)qualities[region][i][j] / 10.0);
			//uint32_t prior3 = ((float) pow((double) 10, -(double)qualities[region][i][j] / 10.0)) << 23;
			
			
			prior_array[i][2 * j] = (int) ((1 - prior) * ONE);
			prior_array[i][2*j + 1] = (int) ((prior / 3)*ONE);
		}
	}
	return prior_array;
}

void print(int** prior, int region) {
	for (int i = 0; i < nr_reads[region]; i++) {
		for (int j = 0; j < reads_len[region][i] - 1; j++) {
			printf("%d %d ", prior[i][2 * j], prior[i][2 * j + 1]);
		}
		printf("\n");
	}
}

void populate_mram(struct dpu_set_t set, uint32_t nr_dpus, int region) {
	
	struct dpu_set_t dpu;
	uint32_t each_dpu;
	uint32_t offset = 0;
	uint32_t offsets[6]; //6 is the number of arrays we need to copy to the mram
	offsets[0] = 0;
	
	
	DPU_ASSERT(dpu_broadcast_to(set, DPU_MRAM_HEAP_POINTER_NAME, 0, reads_len[region], (nr_reads[region] + nr_reads[region] % 2) *sizeof(uint32_t), DPU_XFER_DEFAULT));
	offset += (nr_reads[region] + nr_reads[region] % 2) * sizeof(uint32_t);
	offsets[1] = offset;


	DPU_ASSERT(dpu_broadcast_to(set, DPU_MRAM_HEAP_POINTER_NAME, offset, haplotypes_len[region], (nr_haplotypes[region] + nr_haplotypes[region] % 2) * sizeof(uint32_t), DPU_XFER_DEFAULT));
	offset += (nr_haplotypes[region] + nr_haplotypes[region] % 2) * sizeof(uint32_t);
	offsets[2] = offset;
	

	int reads_arr_size = 0;
	
	/*
	reads_arr_size = total_read_length[region];
	DPU_ASSERT(dpu_broadcast_to(set, DPU_MRAM_HEAP_POINTER_NAME, offset, reads[region], total_read_length[region] * sizeof(char), DPU_XFER_DEFAULT));
	offset += reads_arr_size * sizeof(char);
	*/
	for (int i = 0; i < nr_reads[region]; i++) {
		
		reads_arr_size = (reads_len[region][i] % 8 == 0) ? reads_len[region][i] : reads_len[region][i] + 8 - reads_len[region][i] % 8;
		DPU_ASSERT(dpu_broadcast_to(set, DPU_MRAM_HEAP_POINTER_NAME, offset, reads_array[region][i], reads_arr_size * sizeof(char), DPU_XFER_DEFAULT));
		offset += reads_arr_size * sizeof(char);
	}
	offsets[3] = offset;

	int** prior_array = create_prior(region);
	//print(prior_array, region);
	int prior_arr_size = 0;

	for (int i = 0; i < nr_reads[region]; i++) {
		prior_arr_size = 2 * (reads_len[region][i]) * sizeof(int);
		DPU_ASSERT(dpu_broadcast_to(set, DPU_MRAM_HEAP_POINTER_NAME, offset, prior_array[i], prior_arr_size, DPU_XFER_DEFAULT));
		offset += prior_arr_size;
	}
	offsets[4] = offset;
	free_prior(prior_array, region);

	int hap_arr_size = 0;
	

	/*
	hap_arr_size = total_hap_length[region];
	DPU_ASSERT(dpu_broadcast_to(set, DPU_MRAM_HEAP_POINTER_NAME, offset, haplotypes[region], total_hap_length[region] * sizeof(char), DPU_XFER_DEFAULT));
	offset += hap_arr_size * sizeof(char);
	*/

	for (int i = 0; i < nr_haplotypes[region]; i++) {
		hap_arr_size = (haplotypes_len[region][i] % 8 == 0) ? haplotypes_len[region][i] : haplotypes_len[region][i] + 8 - haplotypes_len[region][i] % 8;
		DPU_ASSERT(dpu_broadcast_to(set, DPU_MRAM_HEAP_POINTER_NAME, offset, haplotypes_array[region][i], hap_arr_size *sizeof(char), DPU_XFER_DEFAULT));
		offset += hap_arr_size * sizeof(char);
	}
	offsets[5] = offset;

	
	DPU_ASSERT(dpu_broadcast_to(set, DPU_MRAM_HEAP_POINTER_NAME, offset, haplotypes_val[region], (nr_haplotypes[region] + nr_haplotypes[region] % 2) * sizeof(uint32_t), DPU_XFER_DEFAULT));
	offset += nr_haplotypes[region] * sizeof(uint32_t);


	DPU_ASSERT(dpu_broadcast_to(set, "heap_offsets", 0, &offsets, sizeof(offsets), DPU_XFER_DEFAULT));
	
	return;
}
