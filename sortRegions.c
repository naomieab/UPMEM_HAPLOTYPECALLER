#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

//This file must contain the array sorted sorted[i]=index of the ith region in increasing order (in term of computation need)
// and defines needed
#include "sortRegions.h"

#define MAX_HAP_NB 4
#define MAX_READ_NB 20//165
#define MAX_COMPLEXITY 5000000

int lines_nb[NR_REGIONS];
char** file[NR_REGIONS];
double** results[NR_REGIONS];
char buffer[BUFFER_SIZE];

int sorted[NR_REGIONS];
int complexity[NR_REGIONS];
int main(int argc, char* argv[]) {
	for (int i = 0; i < NR_REGIONS; i++) { sorted[i] = i; }
	for (int i = 0; i < NR_REGIONS; i++) { complexity[i] = 1; }//removes the complexity condition
	if (argc != 6) {
		printf("Wrong parameters! You must provide 1st argument file to read, 2nd argument file to write regions, 3rd argument file of results, 4th argument file to write new splitted results and 5th argument mode (split/sort)\n");
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

	FILE* resultfile = fopen(argv[3], "r");
	if (!resultfile) {
		printf("Cannot open result file");
		return 0;
	}
	FILE* newresult = fopen(argv[4], "w");
	if (!newresult) {
		printf("Cannot open new result file");
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
			fgets(buffer, BUFFER_SIZE, readfile);
			fgets(buffer, BUFFER_SIZE, readfile);
		}
		file[region] = malloc(lines * sizeof(char*));
		if (!file[region]) {
			return 0;
		}
		results[region] = malloc(reads * sizeof(double*));
		if (!results[region]) { return 0; }
		for (int i = 0; i < reads; i++) {
			results[region][i] = malloc(haps * sizeof(double));
			if (!results[region][i]) { return 0; }
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
		reads = (lines_nb[region] - haps - 2) /2;
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
		//printf("Reads %d and hap %d\n", reads, haps);
		for (i = 0; i < reads; i++) {
			fgets(buffer, BUFFER_SIZE, resultfile);
			//printf("write hap %s\n", buffer);
			char* tmp = strtok(buffer, ",");
			int j = 0;
			while(tmp!=NULL){
				results[region][i][j] = atof(tmp);
				tmp = strtok(NULL, ",");
				j++;
			}
		}
		// 3 blank lines to remove in result file
		fgets(buffer, BUFFER_SIZE, resultfile);
		fgets(buffer, BUFFER_SIZE, resultfile);
		fgets(buffer, BUFFER_SIZE, resultfile);
		region++;
	}


	FILE* newfile = fopen(argv[2], "w");
	int split;
	//long int ref = 7000000;
	counter = 0;
	if (strcmp(argv[3], "sort")==0) {
		printf("Sort mode!\n");
		for (i = 0; i < NR_REGIONS; i++) {
			//printf("REgion %d complexity %d sorted[i] %d\n", i, complexity[sorted[i]], sorted[i]);
			//if ((i + counter) % 2546 == 0) { ref = complexity[sorted[i]]; }
			if (complexity[sorted[i]] < MAX_COMPLEXITY) {
				for (int j = 0; j < lines_nb[sorted[i]]; j++) {
					fprintf(newfile, "%s", file[sorted[i]][j]);
				}
				nb_hap = atoi(file[sorted[i]][0]);
				nb_read = atoi(file[sorted[i]][nb_hap + 1]);
				for (int j = 0; j < nb_read; j++) {
					for (int k = 0; k < nb_hap; k++) {
						fprintf(newresult, "%f,", results[sorted[i]][j][k]);
					}
					fprintf(newresult, "\n");
				}
				fprintf(newresult, "\n\n\n");
			}
			else {
				split = round((double)complexity[sorted[i]] / MAX_COMPLEXITY);
				//printf("COmplexity is %ld and splitting in %d\n", complexity[sorted[i]], split);

				counter += (split - 1);
				//if ((i + counter) % 2546 == 0 || ref ==0) { ref = complexity[sorted[i]]/split; }
				nb_hap = atoi(file[sorted[i]][0]);
				nb_read = atoi(file[sorted[i]][nb_hap + 1]);
				//printf("Ref is %d Split into %d with nb read %d\n", ref, split, nb_read);
				limit = nb_read / split;
				for (int j = 0; j < split; j++) {
					reads = (j == split - 1) ? limit : limit + nb_read % split;
					for (k = 0; k < nb_hap+1; k++) {
						fprintf(newfile, "%s", file[sorted[i]][k]);
					}
					fprintf(newfile, "%d\n", reads);
					for (k = 0; k < reads; k++) {
						fprintf(newfile, "%s", file[sorted[i]][nb_hap + 2 + j*limit +k]);
					}
					for (int l = 0; l < reads; l++) {
						for (int k = 0; k < nb_hap; k++) {
							fprintf(newresult, "%f,", results[sorted[i]][j*limit + l][k]);
						}
						fprintf(newresult, "\n");
					}
					fprintf(newresult, "\n\n\n");
				}
			
			}

		}
		fclose(newfile);
		printf("Number of added regions: %d\n", counter);
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
						//if (read_limit > 1) { printf("The string is: %s", file[sorted[i]][nb_hap + 2 + 2 * l * MAX_READ_NB + r1]); }
						fprintf(newfile, "%s", file[sorted[i]][nb_hap + 2 + 2*l * MAX_READ_NB + r1]);
						
					}
					for (int r1 = 0; r1 < r; r1++) {
						for (int h1 = 0; h1 < h; h1++) {
							fprintf(newresult, "%f,", results[sorted[i]][l * MAX_READ_NB + r1 ][k * MAX_HAP_NB + h1]);
						}
						fprintf(newresult, "\n");
					}
					//fprintf(newresult, "\n\n\n");
				}
			}
			counter = counter + (limit * read_limit) - 1;
	}
	printf("Number of added regions:%d\n", counter);
	return 1;
}



