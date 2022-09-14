#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

//This file must contain the array sorted sorted[i]=index of the ith region in increasing order (in term of computation need)
// and defines needed
#include "sortRegions.h"

#define MAX_HAP_NB 24
#define MAX_READ_NB 120//200

int lines_nb[NR_REGIONS];
char** file[NR_REGIONS];
char buffer[BUFFER_SIZE];

int sorted[NR_REGIONS];

int main(int argc, char* argv[]) {
	for (int i = 0; i < NR_REGIONS; i++) { sorted[i] = i; }
	if (argc != 4) {
		printf("Wrong parameters! You must provide 1st argument file to read and 2nd argument file to write regions and 3rd argument mode (split/sort)\n");
		return 0;
	}

	FILE* readfile = fopen(argv[1], "r");
	if (!readfile) {
		printf("Cannot open input file");
		return 0;
	}

	FILE* sortedfile = fopen(argv[2], "w");
	if (!sortedfile) {
		printf("Cannot open output file");
		return 0;
	}
	
	int haps, reads, lines, region = 0;
	int i, len, counter=0;
	char* str;
	int nb_hap, k, limit;
	int nb_read, l, read_limit;
	while (fgets(buffer, BUFFER_SIZE, readfile) != 0) {
		haps = atoi(buffer);
		assert(haps > 0);
		lines = haps + 1;
		for (i = 0; i < haps; i++) {
			fgets(buffer, BUFFER_SIZE, readfile);
		}
		reads = atoi(fgets(buffer, BUFFER_SIZE, readfile));
		assert(reads > 0);
		lines += (2*reads + 1);
		for (i = 0; i < reads; i++) {
			fgets(buffer, BUFFER_SIZE, readfile); //get read line
			fgets(buffer, BUFFER_SIZE, readfile); //get transitions line
		}
		file[region] = malloc(lines * sizeof(char*));
		if (!file[region]) {
			return 0;
		}
		lines_nb[region] = lines;
		region++;
	}
	fclose(readfile);
	
	readfile = fopen(argv[1], "r");
	if (!readfile) {
		printf("Cannot open input file");
		return 0;
	}

	region = 0;
	while (fgets(buffer, BUFFER_SIZE, readfile) != 0 && region < NR_REGIONS) {
		str = buffer;
		file[region][0] = malloc((strlen(str) + 1) * sizeof(char));
		if (!file[region][0]) { printf("ERROR!\n"); return 0; }
		haps = atoi(buffer);
		strcpy(file[region][0], str);
		assert(atoi(buffer) > 0);
		for (i = 1; i < haps + 1; i++) {
			fgets(buffer, BUFFER_SIZE, readfile);
			str = buffer;
			file[region][i] = malloc((strlen(str) + 1) * sizeof(char));
			if (!file[region][i]) { printf("ERROR!\n"); return 0; }

			strcpy(file[region][i], str);
		}
		reads = (lines_nb[region] - haps - 2)/2;
		fgets(buffer, BUFFER_SIZE, readfile);
		str = buffer;
		file[region][i] = malloc((strlen(str) + 1) * sizeof(char));
		if (!file[region][i]) { printf("ERROR!\n"); return 0; }
		strcpy(file[region][i], str);
		i++;
		for (i; i < haps + 2*reads + 2; i++) {
			fgets(buffer, BUFFER_SIZE, readfile);
			str = buffer;
			file[region][i] = malloc((strlen(str) + 1) * sizeof(char));
			if (!file[region][i]) { printf("ERROR!\n"); return 0; }
			strcpy(file[region][i], str);
		}
		region++;
	}


	FILE* newfile = fopen(argv[2], "w");

	if (strcmp(argv[3], "sort")==0) {
		printf("Sort mode!\n");
		for (i = 0; i < NR_REGIONS; i++) {
			for (int j = 0; j < lines_nb[sorted[i]]; j++) {
				fprintf(newfile, "%s", file[sorted[i]][j]);
			}
		}
		fclose(newfile);
		return 1;
	}



	for (i = 0; i < NR_REGIONS; i++) {
		assert(i < NR_REGIONS); 
		assert(sorted[i] < NR_REGIONS);
		nb_hap = atoi(file[sorted[i]][0]);
		assert(nb_hap != 0);
			limit = ceil((double)nb_hap / MAX_HAP_NB);
			nb_read = atoi(file[sorted[i]][nb_hap + 1]);
			assert(nb_read != 0);
				read_limit = ceil((double)nb_read / MAX_READ_NB);
			for ( k = 0; k < limit; k++) {
				int h = (k == limit-1) ? nb_hap % MAX_HAP_NB : MAX_HAP_NB;
				if (h == 0) { h = MAX_HAP_NB; } //case when the number of haps is a multiple of MAX_HAP_NB
				for (l = 0; l < read_limit; l++) {
					int r = (l == read_limit - 1) ? nb_read % MAX_READ_NB : MAX_READ_NB;
					if (r == 0) { r = MAX_READ_NB; } //case when the number of haps is a multiple of MAX_READ_NB
					fprintf(newfile, "%d\n", h);
					//write haplotype sequence
					assert(i < NR_REGIONS);
					assert(sorted[i] < NR_REGIONS);

					for (int h1 = 0; h1 < h; h1++) {
						fprintf(newfile, "%s", file[sorted[i]][k * MAX_HAP_NB + 1 + h1]);
					}
					assert(i < NR_REGIONS);
					fprintf(newfile, "%d\n", r);
					for (int r1 = 0; r1 < 2*r; r1++) {
						fprintf(newfile, "%s", file[sorted[i]][nb_hap + 2 + l * 2 * MAX_READ_NB + r1]);
					}
				}
			}
			counter = counter + (limit * read_limit) - 1;
	}
	printf("Number of added regions:%d\n", counter);
	return 1;
}



