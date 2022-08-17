#include "parser.h"
#include "constants.h"


uint32_t offset[NR_REGIONS][OFFSET_SIZE];


uint32_t nr_haplotypes[NR_REGIONS]; //an array keeping number of haplotypes in all regions
uint32_t nr_reads[NR_REGIONS]; //idem as haplotypes


//TOTAL_READS= maximum total number of reads in a chunk of NR_DPUS regions + MAX_NUMBER OF READS IN ONE REGION  
uint32_t reads_len[TOTAL_READS];
//TOTAL_READS_SIZE = TOTAL_READS * MAX READ LEN  
char reads_array[TOTAL_READS * MAX_READ_LENGTH];
uint32_t haplotypes_len[TOTAL_HAPS];
uint32_t haplotypes_val[TOTAL_HAPS];
char haplotypes_array[TOTAL_HAPS * MAX_HAPLOTYPE_LENGTH];
uint32_t priors[TOTAL_READS * MAX_READ_LENGTH * 2];



//hap_idx is the partial sum of the all the reads in the regions before
//hap_idx+index is the offset where to write this current haplotype
void add_haplotype(FILE* file, int hap_idx, int index) {
	if (hap_idx + index > TOTAL_HAPS) { fprintf(stderr, "Error, number of haps %d is bigger than allocated with TOTAL_HAPS\n", hap_idx + index); }
	char buffer[BUFFER_SIZE];
	assert(fgets(buffer, BUFFER_SIZE, file));
	char* hap_str = strtok(buffer, ",");
	int hap_length = strlen(hap_str);
	haplotypes_val[hap_idx + index] = (int)(log10(1.0 / (hap_length)) * ONE);
	haplotypes_len[hap_idx + index] = hap_length;
	strncpy(&haplotypes_array[(hap_idx + index) * MAX_HAPLOTYPE_LENGTH], hap_str, MAX_HAPLOTYPE_LENGTH);
	return;
}



//read_idx is the partial sum of the all the reads in the regions before
void add_read(FILE* file, int read_idx, int index) {
	if (read_idx + index > TOTAL_READS) { fprintf(stderr, "Error, number of reads %d is bigger than allocated with TOTAL_READS\n", read_idx + index); }
	char buffer[BUFFER_SIZE];
	assert(fgets(buffer, BUFFER_SIZE, file));
	char* token = strtok(buffer, ",");
	int read_length = strlen(token);
	reads_len[read_idx + index] = read_length;
	strncpy(&reads_array[(read_idx + index) * MAX_READ_LENGTH], token, MAX_READ_LENGTH);
	int j = 0;
	for (token = strtok(NULL, ","); token != NULL && j < read_length; token = strtok(NULL, ","), j++) {
		int quality = atoi(token);
		double probLog10 = log10(1 - pow((double)10, -(double)quality / 10.0));
		double errorProbLog10 = log10(pow((double)10, -(double)quality / 10.0)) - log10(3);
		priors[(read_idx + index) * 2 * MAX_READ_LENGTH + 2 * j] = (int)(probLog10 * ONE);
		priors[(read_idx + index) * 2 * MAX_READ_LENGTH + 2 * j + 1] = (int)(errorProbLog10 * ONE);
	}
}




FILE* read_data(FILE* file, int nr_dpus) {
	//FILE* file = fopen(filename, "r");
	if (!file) {
		printf("Wrong input file");
		return NULL;
	}
	int current_region = 0;
	int hap_idx = 0, read_idx = 0;
	char buffer[BUFFER_SIZE];
	while (current_region < nr_dpus && fgets(buffer, BUFFER_SIZE, file) != 0) {
		nr_haplotypes[current_region] = atoi(buffer);
		for (int i = 0; i < nr_haplotypes[current_region]; i++) {
			add_haplotype(file, hap_idx, i);
		}
		hap_idx += nr_haplotypes[current_region];
		nr_reads[current_region] = atoi(fgets(buffer, BUFFER_SIZE, file));
		for (int i = 0; i < nr_reads[current_region]; i++) {
			add_read(file, read_idx, i);
		}
		read_idx += nr_reads[current_region];
		current_region++;
	}
	for (int i = 1; i < NR_REGIONS; i++) {
		offset[i][READS_LEN_ARRAY] = offset[i - 1][READS_LEN_ARRAY] + nr_reads[i - 1];// + (nr_reads[i-1]%2==1);
		offset[i][HAPLOTYPES_LEN_VAL_ARRAY] = offset[i - 1][HAPLOTYPES_LEN_VAL_ARRAY] + nr_haplotypes[i - 1]; //+ (nr_haplotypes[i-1]%2==1);
		offset[i][READS_ARR] = offset[i - 1][READS_ARR] + (nr_reads[i - 1] * MAX_READ_LENGTH);
		offset[i][HAPS_ARR] = offset[i - 1][HAPS_ARR] + (nr_haplotypes[i - 1] * MAX_HAPLOTYPE_LENGTH);
		offset[i][PRIOR_ARR] = offset[i - 1][PRIOR_ARR] + (nr_reads[i - 1] * MAX_READ_LENGTH * 2);
	}
	return file;
}




void free_mem(FILE* file) {
	fclose(file);
}
