#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define BUFFER_SIZE 2000
#define NR_REGIONS 132108

#define MAX_HAP_NB 4
#define MAX_READ_NB 20//165

#define FILE_MAX_READS 6881
#define FILE_MAX_HAPLOTYPES 165
#define FILE_MAX_READ_LEN 80 //must be the maximum read length + 1 (because the line contains 1 str + read len qualities)
#define FILE_MAX_HAP_LEN 373


char buffer[BUFFER_SIZE];

char reads_buffer[2*FILE_MAX_READS][FILE_MAX_READ_LEN];
char haplotypes[FILE_MAX_HAPLOTYPES][FILE_MAX_HAP_LEN];
double results[FILE_MAX_READS][FILE_MAX_HAPLOTYPES];

int hap_split[NR_REGIONS];
int read_split[NR_REGIONS];
int hap_number[NR_REGIONS];
int read_number[NR_REGIONS];

int main(int argc, char* argv[]) {
	if (argc != 6) {
		printf("Wrong parameters! You must provide 1st argument file to read, 2nd argument file to write regions, 3rd argument file of results, 4th argument file to write new splitted results and 5th argument mode (split/sort)\n");
		return 0;
	}

	FILE* readfile = fopen(argv[1], "r");
	if (!readfile) {
		printf("Cannot open input file");
		return 0;
	}

	FILE* newfile = fopen(argv[2], "w");
	if (!newfile) {
		printf("Cannot open output file");
		return 0;
	}

	FILE* resultfile = fopen(argv[3], "r");
	if (!resultfile) {
		printf("Cannot open result file");
		return 0;
	}
	FILE* newresult = fopen(argv[4], "w");
	if (!newresult) {
		printf("Cannot open output result file");
		return 0;
	}

	int haps, reads, region = 0;
	int i, len, counter = 0;
	char* str;
	int k, limit;
	int l, read_limit;


	while (fgets(buffer, BUFFER_SIZE, readfile) != 0 && region < NR_REGIONS) {
		//Read current input region
		haps = atoi(buffer);
		assert(haps > 0);
		for (i = 0; i < haps; i++) {
			fgets(buffer, BUFFER_SIZE, readfile); //read current haplotype string
			strcpy(haplotypes[i], buffer);
		}
		fgets(buffer, BUFFER_SIZE, readfile); //read number of reads
		reads = atoi(buffer);
		for (i = 0; i < reads; i++) {
			fgets(buffer, BUFFER_SIZE, readfile); //read the current read string + qualities
			strcpy(reads_buffer[2 * i], buffer);
			fgets(buffer, BUFFER_SIZE, readfile); //read the transition qualities
			strcpy(reads_buffer[2 * i + 1], buffer);
		}
		//Read current region results
		for (i = 0; i < reads; i++) {
			assert(fgets(buffer, BUFFER_SIZE, resultfile));
			char* tmp = strtok(buffer, ",");
			int j = 0;
			while (j < haps) {
				assert(tmp != NULL);
				results[i][j] = atof(tmp);
				tmp = strtok(NULL, ",");
				j++;
			}
		}
		// 3 blank lines to remove in result file
		assert(fgets(buffer, BUFFER_SIZE, resultfile));
		assert(fgets(buffer, BUFFER_SIZE, resultfile));
		assert(fgets(buffer, BUFFER_SIZE, resultfile));
		region++;
		printf("Region %d\n",region);
		//Write split region to file
		limit = ceil((double)haps / MAX_HAP_NB);
		read_limit = ceil((double)reads / MAX_READ_NB);

		//prepare data for post-processing
		hap_split[region] = limit;
		read_split[region] = read_limit;
		hap_number[region] = haps;
		read_number[region] = reads;

		for (l = 0; l < read_limit; l++) {
			int r = (l == read_limit - 1) ? reads % MAX_READ_NB : MAX_READ_NB;
			if (r == 0) { r = MAX_READ_NB; } //Cover case when the number of haps is a multiple of MAX_READ_NB

			for (k = 0; k < limit; k++) {
				int h = (k == limit - 1) ? haps % MAX_HAP_NB : MAX_HAP_NB;
				if (h == 0) { h = MAX_HAP_NB; } //Cover case when the number of haps is a multiple of MAX_HAP_NB
				fprintf(newfile, "%d\n", h); //Write number of haplotypes
				for (int h1 = 0; h1 < h; h1++) {
					fprintf(newfile, "%s", haplotypes[k * MAX_HAP_NB + h1]); //Write haplotype sequence
				}
				fprintf(newfile, "%d\n", r); //Write number of reads
				for (int r1 = 0; r1 < 2 * r; r1++) {
					fprintf(newfile, "%s", reads_buffer[2 * l * MAX_READ_NB + r1]); //Write read/transition sequence
				}
				//Write result to result file
				for (int r1 = 0; r1 < r; r1++) {
					for (int h1 = 0; h1 < h; h1++) {
						fprintf(newresult, "%f,", results[l * MAX_READ_NB + r1][k * MAX_HAP_NB + h1]);
					}
					fprintf(newresult, "\n");
				}
			}
		}
		counter = counter + (limit * read_limit) - 1;
	}
	fclose(readfile);
	fclose(resultfile);
	printf("Number of added regions:%d\n", counter);
	return 1;
	post_processing();
}

void post_processing(FILE* split_output, FILE* output) {
	int region = 0;
	char* tmp;
	while (region < NR_REGIONS) {
		for (int i = 0; i < read_split[region]; i++) {
			for (int j = 0; j < hap_split[region]; j++) {
				int reads = (i == (read_split[region] - 1)) : read_number[region] % MAX_READ_NB : MAX_READ_NB;
				int haps = (j == (hap_split[region] - 1)) : hap_number[region] % MAX_HAP_NB : MAX_HAP_NB;
				for (int k = 0; k < reads; k++) {
					fgets(buffer, BUFFER_SIZE, split_output);
					tmp = strtok(buffer, ",");
					for (int l = 0; l < haps; l++) {
						results[i * MAX_READ_NB + k][j * MAX_HAP_NB + l] = atof(tmp);
						tmp = strtok(NULL, ",");
					}
				}

			}
		}
		for (int i = 0; i < read_number[region]; i++) {
			for (int j = 0; j < hap_number[region]; j++) {
				fprintf(output, "%s,", results[i][j]);
			}
			fprintf(output, "\n");
		}
		fprintf(output, "\n\n\n");
	}
}