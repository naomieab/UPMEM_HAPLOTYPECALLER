#include "parser.h"
#include "constants.h"

uint32_t nr_regions; //number of regions


/*
uint32_t** reads_len;
char*** reads_array;
uint32_t*** qualities;
uint32_t** haplotypes_len;
uint32_t** haplotypes_val;
char*** haplotypes_array;
*/
char** reads;
char** haplotypes;

uint32_t nr_haplotypes[NR_REGIONS]; //an array keeping number of haplotypes in all regions
uint32_t nr_reads[NR_REGIONS]; //idem as haplotypes

uint32_t reads_len[NR_REGIONS][MAX_READ_NUM];
char reads_array[NR_REGIONS][MAX_READ_NUM][MAX_READ_LENGTH];

uint32_t qualities[NR_REGIONS][MAX_READ_NUM][2*MAX_READ_LENGTH];
uint32_t haplotypes_len[NR_REGIONS][MAX_HAPLOTYPE_NUM];
uint32_t haplotypes_val[NR_REGIONS][MAX_HAPLOTYPE_NUM];
char haplotypes_array[NR_REGIONS][MAX_HAPLOTYPE_NUM][MAX_HAPLOTYPE_LENGTH];



void read_and_allocate(FILE* file) {
	char buffer[BUFFER_SIZE];
	assert(fgets(buffer, BUFFER_SIZE, file) != 0);
	
	nr_regions = atoi(buffer);
	//EVERYTHING IS MADE CONSTANT
	/*
	nr_haplotypes = malloc(nr_regions * sizeof(uint32_t));
	nr_reads = malloc(nr_regions * sizeof(uint32_t));
	reads_len = malloc(nr_regions * sizeof(uint32_t*));
	reads_array = malloc(nr_regions * sizeof(char**));
	reads = malloc(nr_regions * sizeof(char*));
	qualities = malloc(nr_regions * sizeof(uint32_t**));
	haplotypes_len = malloc(nr_regions * sizeof(uint32_t*));
	haplotypes_val = malloc(nr_regions * sizeof(uint32_t*));
	haplotypes_array = malloc(nr_regions * sizeof(char**));
	haplotypes = malloc(nr_regions * sizeof(char*));
	total_read_length = malloc(nr_regions * sizeof(uint32_t));
	total_hap_length = malloc(nr_regions * sizeof(uint32_t));

	assert(nr_haplotypes);
	assert(nr_reads);
	assert(reads_len);
	assert(reads_array);
	assert(reads);
	assert(qualities);
	assert(haplotypes_len);
	assert(haplotypes_val);
	assert(haplotypes_array);
	assert(haplotypes);
	assert(total_read_length);
	assert(total_hap_length);
	*/
}


void add_haplotype(FILE* file, int region, int index) {
	char buffer[BUFFER_SIZE];
	assert(fgets(buffer, BUFFER_SIZE, file));
	char* hap_str = strtok(buffer, ",");
	int hap_length = strlen(hap_str);
	haplotypes_len[region][index] = hap_length;
	haplotypes_val[region][index] = (int) (log(1.0 / (hap_length + 1)) * ONE);
	//could erase the following lines since we work now with constant numbers
	/*total_hap_length[region] += hap_length;
	if (hap_length % 8 != 0) {
		hap_length += 8 - hap_length % 8;
	}
	haplotypes_array[region][index] = malloc(MAX_HAPLOTYPE_LENGTH  * sizeof(char));
	*/
	strncpy(haplotypes_array[region][index], hap_str, MAX_HAPLOTYPE_LENGTH/*hap_length*/);
}

void add_read(FILE* file, int region, int index) {
	char buffer[BUFFER_SIZE];
	assert(fgets(buffer, BUFFER_SIZE, file));
	char* token = strtok(buffer, ",");
	int read_length = strlen(token);
	reads_len[region][index] = read_length;
	//could erase the following lines since we work now with constant numbers
	/*int malloc_read_length = read_length;
	if (malloc_read_length % 8 != 0) {
		malloc_read_length += 8 - malloc_read_length % 8;
	}
	total_read_length[region] += read_length;
	reads_array[region][index] = malloc(MAX_READ_LENGTH * sizeof(char));
	qualities[region][index] = malloc(MAX_READ_LENGTH * sizeof(uint32_t));*/
	strncpy(reads_array[region][index], token, MAX_READ_LENGTH/*read_length*/);
	int j = 0;
	for(token=strtok(NULL,","); token != NULL && j < read_length; token = strtok(NULL, ","), j++) {
		qualities[region][index][j] = atoi(token);
	}	
}

/*
void concatenate_reads(int region) {
	if (total_read_length[region] % 8 != 0) {
		total_read_length[region] += 8 - (total_read_length[region] % 8);
	}
	reads[region] = malloc(total_read_length[region] * sizeof(char));
	assert(reads[region]);
	int tmp = 0;
	for (int i = 0; i < nr_reads[region]; i++) {
		strncpy(&reads[region][tmp], reads_array[region][i], reads_len[region][i]);
		tmp += reads_len[region][i];
	}
}

void concatenate_haplotypes(int region) {
	if (total_hap_length[region] % 8 != 0) {
		total_hap_length[region] += (8 - total_hap_length[region] % 8);
	}
	haplotypes[region] = malloc(total_hap_length[region] * sizeof(char));
	assert(haplotypes[region]);
	int tmp = 0;
	for (int i = 0; i < nr_haplotypes[region]; i++) {

		strncpy(&haplotypes[region][tmp], haplotypes_array[region][i], haplotypes_len[region][i]);
		tmp += haplotypes_len[region][i];
	}
}
*/

FILE* read_data(char* filename) {
	FILE* file = fopen(filename, "r");
	if (!file) {
		printf("Oupsi! Can't read input file");
		return NULL;
	}
	read_and_allocate(file);
	int current_region = 0;
	char buffer[BUFFER_SIZE];
	while (fgets(buffer, BUFFER_SIZE, file) != 0 && current_region < NR_REGIONS) {
		nr_haplotypes[current_region] = atoi(buffer);

		//erase following line
		//nr_haplotypes[current_region] = nr_haps;
		//haplotypes_len[current_region] = malloc( MAX_HAPLOTYPE_NUM * sizeof(uint32_t));
		//haplotypes_val[current_region] = malloc( MAX_HAPLOTYPE_NUM * sizeof(uint32_t));
		//haplotypes_array[current_region] = malloc(MAX_HAPLOTYPE_NUM * sizeof(char*));
		//total_hap_length[current_region] = 0;
		for (int i = 0; i < nr_haplotypes[current_region]; i++) {
			add_haplotype(file, current_region, i);
		}
		//printf("\n HAP LENGTH=%d\n", )
		nr_reads[current_region] = atoi(fgets(buffer, BUFFER_SIZE, file));
		//erase following line
		//nr_reads[current_region] = nr_read;
		//reads_len[current_region] = malloc( MAX_READ_NUM * sizeof(uint32_t));
		//reads_array[current_region] = malloc( MAX_READ_NUM * sizeof(char*));
		//qualities[current_region] = malloc( MAX_READ_NUM * sizeof(uint32_t*));
		//total_read_length[current_region] = 0;
		for (int i = 0; i < nr_reads[current_region]; i++) {
			add_read(file, current_region, i);
		}

		//TODO: change the code to create directly an array containing all the reads together
		//Currently doing copy paste below
		//Idem for haplotypes
		
		//concatenate_reads(current_region);
		//concatenate_haplotypes(current_region);

		//nr_haplotypes[current_region] = nr_haps;
		//nr_reads[current_region] = nr_read;
		printf("Region %d allocated: nr_haps = %d and nr_read = %d\n", current_region, nr_haplotypes[current_region], nr_reads[current_region]);
		current_region++;
	}


	/*for (int i = 0; i < NR_REGIONS; i++) {
		for (int j = 0; j < MAX_HAPLOTYPE_NUM; j++) {
			printf("%s\n", haplotypes_array[i][j]);
		}
		printf("\n\n");
	}*/
	return file;
}

void free_mem(FILE* file) {
	int nr_read, nr_haps;
	/*
	for (int i = 0; i < nr_regions; i++) {
		nr_read = nr_reads[i];
		nr_haps = nr_haplotypes[i];
		for (int j = 0; j < nr_read; j++) {
			free(reads_array[i][j]);
			free(qualities[i][j]);
		}
		for (int k = 0; k < nr_haps; k++) {
			free(haplotypes_array[i][k]);
		}
		free(reads_array[i]);
		free(reads[i]);
		free(qualities[i]);
		free(reads_len[i]);
		free(haplotypes_len[i]);
		free(haplotypes_val[i]);
		free(haplotypes_array[i]);
		free(haplotypes[i]);
	}
	free(reads_array);
	free(reads);
	free(qualities);
	free(reads_len);
	free(haplotypes_len);
	free(haplotypes_array);
	free(haplotypes_val);
	free(haplotypes);
	
	free(total_read_length);
	free(total_hap_length);

	free(nr_haplotypes);
	free(nr_reads);
	*/
	fclose(file);

}