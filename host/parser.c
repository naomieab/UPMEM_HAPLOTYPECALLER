#include "parser.h"
#include "constants.h"

uint32_t nr_regions; //number of regions




uint32_t nr_haplotypes[NR_REGIONS]; //an array keeping number of haplotypes in all regions
uint32_t nr_reads[NR_REGIONS]; //idem as haplotypes

uint32_t reads_len[NR_REGIONS][MAX_READ_NUM];
char reads_array[NR_REGIONS][MAX_READ_NUM][MAX_READ_LENGTH];
uint32_t qualities[NR_REGIONS][MAX_READ_NUM][2*MAX_READ_LENGTH];
uint32_t haplotypes_len[NR_REGIONS][MAX_HAPLOTYPE_NUM];
uint32_t haplotypes_val[NR_REGIONS][MAX_HAPLOTYPE_NUM];
char haplotypes_array[NR_REGIONS][MAX_HAPLOTYPE_NUM][MAX_HAPLOTYPE_LENGTH];



void add_haplotype(FILE* file, int region, int index) {
	char buffer[BUFFER_SIZE];
	assert(fgets(buffer, BUFFER_SIZE, file));
	char* hap_str = strtok(buffer, ",");
	int hap_length = strlen(hap_str);
	haplotypes_len[region][index] = hap_length;
	haplotypes_val[region][index] = (int) (log(1.0 / (hap_length + 1)) * ONE);
	strncpy(haplotypes_array[region][index], hap_str, MAX_HAPLOTYPE_LENGTH/*hap_length*/);
}




void add_read(FILE* file, int region, int index) {
	char buffer[BUFFER_SIZE];
	assert(fgets(buffer, BUFFER_SIZE, file));
	char* token = strtok(buffer, ",");
	int read_length = strlen(token);
	reads_len[region][index] = read_length;
	strncpy(reads_array[region][index], token, MAX_READ_LENGTH);
	int j = 0;
	for(token=strtok(NULL,","); token != NULL && j < read_length; token = strtok(NULL, ","), j++) {
		qualities[region][index][j] = atoi(token);
	}	
}




FILE* read_data(char* filename) {
	FILE* file = fopen(filename, "r");
	if (!file) {
		printf("Oupsi! Can't read input file");
		return NULL;
	}
	int current_region = 0;
	char buffer[BUFFER_SIZE];
	while (fgets(buffer, BUFFER_SIZE, file) != 0 && current_region < NR_REGIONS) {
		nr_haplotypes[current_region] = atoi(buffer);
		for (int i = 0; i < nr_haplotypes[current_region]; i++) {
			add_haplotype(file, current_region, i);
		}
		nr_reads[current_region] = atoi(fgets(buffer, BUFFER_SIZE, file));
		for (int i = 0; i < nr_reads[current_region]; i++) {
			add_read(file, current_region, i);
		}
		printf("Region %d allocated: nr_haps = %d and nr_read = %d\n", current_region, nr_haplotypes[current_region], nr_reads[current_region]);
		current_region++;
	}
	return file;
}




void free_mem(FILE* file) {
	fclose(file);
}