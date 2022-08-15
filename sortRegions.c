#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>


//This file must contain the array sorted sorted[i]=index of the ith region in increasing order (in term of computation need)
// and defines needed
#include "sortRegions.h"


int lines_nb[NR_REGIONS];
char** file[NR_REGIONS];
char buffer[BUFFER_SIZE];


int main(int argc, char* argv[]) {
	if (argc != 3) {
		printf("Wrong parameters! You must provide 1st argument file to read and 2nd argument file to write regions\n");
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
	int i;
	while (fgets(buffer, BUFFER_SIZE, readfile) != 0) {
		haps = atoi(buffer);
		lines = haps + 1;
		for (i = 0; i < haps; i++) {
			fgets(buffer, BUFFER_SIZE, readfile);
		}
		reads = atoi(fgets(buffer, BUFFER_SIZE, readfile));
		lines += (reads + 1);
		for (i = 0; i < reads; i++) {
			fgets(buffer, BUFFER_SIZE, readfile);
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
	int len;
	char* str;
	while (fgets(buffer, BUFFER_SIZE, readfile) != 0 && region < NR_REGIONS) {
		str = buffer;
		file[region][0] = malloc((strlen(str) + 1) * sizeof(char));
		if (!file[region][0]) { printf("ERROR!\n"); return 0; }

		strcpy(file[region][0], str);
		for (i = 1; i < haps + 1; i++) {
			fgets(buffer, BUFFER_SIZE, readfile);
			str = buffer;
			file[region][i] = malloc((strlen(str) + 1) * sizeof(char));
			if (!file[region][i]) { printf("ERROR!\n"); return 0; }

			strcpy(file[region][i], str);
		}
		reads = lines_nb[region] - haps - 2;
		fgets(buffer, BUFFER_SIZE, readfile);
		str = buffer;
		file[region][i] = malloc((strlen(str) + 1) * sizeof(char));
		if (!file[region][i]) { printf("ERROR!\n"); return 0; }
		strcpy(file[region][i], str);
		i++;
		for (i; i < haps + reads + 2; i++) {
			fgets(buffer, BUFFER_SIZE, readfile);
			str = buffer;
			file[region][i] = malloc((strlen(str) + 1) * sizeof(char));
			if (!file[region][i]) { printf("ERROR!\n"); return 0; }

			strcpy(file[region][i], str);
		}
		region++;
	}


	FILE* newfile = fopen(argv[2], "w");
	for (i = 0; i < NR_REGIONS - TO_SEPARATE; i++) {
		for (int j = 0; j < lines_nb[sorted[i]]; j++) {
			fprintf(newfile, "%s", file[sorted[i]][j]);
		}
	}

	/*
	//separation of regions
	//first according to haplotypes
	int counter = 0;
	for (i; i < NR_REGIONS; i++) {
		int haps = atoi(file[sorted[i]][0]);
		for (int j = 0; j < haps; j++) {
			fprintf(newfile, "1\n");
			fprintf(newfile, "%s", file[sorted[i]][j + 1]);
			//don't write first line since it is the number of haplotypes so it is 1 only and also the haplotypes lines
			for (int k = haps + 1; k < lines_nb[sorted[i]]; k++) {
				fprintf(newfile, "%s", file[sorted[i]][k]);
			}
			counter++;
		}
	}*/

	//separation of regions
	//first according to haplotypes 
	//and twice according to reads
	int counter = 0;
	int new_regions_cnt = 0, new_regions_reads;
	int nb_reads, offset, region_reads;

	int max_reads = atoi(file[sorted[i]][haps + 1]); //number of reads in the current region
	for (i; i < NR_REGIONS; i++) {
		if (i % 2546 == 0) {
			max_reads = atoi(file[sorted[i]][haps + 1]);
		}
		int haps = atoi(file[sorted[i]][0]);
		region_reads = atoi(file[sorted[i]][haps + 1]);
		for (int r = 0; r < 2; r++) {
			for (int j = 0; j < haps; j++) {
				if (new_regions_cnt == 0) {
					new_regions_reads = reads;
				}
				//write line with "1" for number of haplotypes
				fprintf(newfile, "1\n");
				//write haplotype sequence 
				fprintf(newfile, "%s", file[sorted[i]][j + 1]);
				//nb_reads is the number of reads to write for the current region
				//offset is the offset where to start the reads (depends on round 0 or round 1)
				//region_reads is the number of reads in the region

				if (r == 0) {
					nb_reads = (int)ceil((double)region_reads / 2);
					offset = 0;
				}
				else if (r == 1) {
					nb_reads = (int)floor((double)region_reads / 2);
					offset = (int)ceil((double)region_reads / 2);
				}
				//write number of reads
				fprintf(newfile, "%d\n", nb_reads);
				//don't write first line since it is the number of haplotypes so it is 1 only and also the haplotypes lines
				for (int k = 0; k < nb_reads; k++) {
					//for (int k = haps + 1; k < lines_nb[sorted[i]]; k++) {
					fprintf(newfile, "%s", file[sorted[i]][k + haps + 2 + offset]);
				}
				counter++;
			}
		}
	}

	printf("Number of added regions:%d\n", counter);
}



