#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <unistd.h> 
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#define BUFFER_SIZE 1024
#define MAX_READ_NUMBER 20000
#define MAX_HAPLOTYPE_NUMBER 500

int add_read(FILE* file);
int add_haplotype(FILE* file);


int reads_len[MAX_READ_NUMBER];
int haplotypes_len[MAX_HAPLOTYPE_NUMBER];

int main(int argc, char* argv[]) {
	if (argc != 3) {
		printf("Wrong parameters! You must provide 1st argument file to read and 2nd argument file to write regions\n");
		return 0;
	}

	FILE* file = fopen(argv[1], "r");
	if (!file) {
		printf("Cannot open input file");
		return 0;
	}

	FILE* stat_file = fopen(argv[2], "w");
	if (!stat_file) {
		printf("Cannot open output file");
		return 0;
	}

	int current_region = 0, nr_haplotypes = 0, nr_reads = 0;
	int reads, haps;
	int max_hap_len = 0, max_hap_nb = 0;
	int max_read_len = 0, max_reads_nb = 0;
	long int complexity_stat = 0;
	int  cnt_read_max = 0, cnt_read=0, counter=0;
	char buffer[BUFFER_SIZE];
	int max_complexity = 0;

	while (fgets(buffer, BUFFER_SIZE, file) != 0 ) {
		haps = atoi(buffer);
		if (haps > max_hap_nb) {
			max_hap_nb = haps;
		}
		nr_haplotypes += (haps + (haps % 2 == 1));
		assert(haps < MAX_HAPLOTYPE_NUMBER);
		for (int i = 0; i < haps; i++) {
			haplotypes_len[i] = add_haplotype(file);
			if (haplotypes_len[i] > max_hap_len) {
				max_hap_len = haplotypes_len[i];
			}
		}
		

		reads = atoi(fgets(buffer, BUFFER_SIZE, file));
		if (reads > max_reads_nb) {
			max_reads_nb = reads;
		}
		nr_reads += (reads + (reads % 2 == 1));

		assert(reads < MAX_READ_NUMBER);
		for (int i = 0; i < reads; i++) {
			reads_len[i] = add_read(file);
			if (reads_len[i] > max_read_len) {
				max_read_len = reads_len[i];
			}
		}

		complexity_stat = 0;
		int avg_hap_len=0;
		for (int i = 0; i < haps; i++) {
			for (int j = 0; j < reads; j++) {
				complexity_stat += haplotypes_len[i] * reads_len[j];
			}
			avg_hap_len += haplotypes_len[i];
		}
		if (complexity_stat > max_complexity) { max_complexity = complexity_stat; }
		fprintf(stat_file,"%d, %ld, %d, %d, %d\n", current_region, complexity_stat, haps, reads, avg_hap_len/haps);
		current_region++;

		counter++;
		cnt_read += reads;
		if (counter % 2546 == 0) {
			if (cnt_read > cnt_read_max) {
				cnt_read_max = cnt_read;
			}
			cnt_read = 0;
			counter = 0;	
		}

	}
	printf("Statistics of the file:\n");
	printf("Number of regions: %d\n", current_region);
	printf("Average number of haplotypes per region: %f\n", (double)nr_haplotypes/current_region);
	printf("Average number of reads per region: %f\n", (double)nr_reads/current_region);
	printf("Max number of haplotypes in region: %d\n", max_hap_nb);
	printf("Max haplotype length: %d\n", max_hap_len);
	printf("Max number of reads in region: %d\n", max_reads_nb);
	printf("Max read length: %d\n", max_read_len);
	printf("Max number of reads in chunck of 2546: %d\n", cnt_read_max);
	printf("Max complexity: %d\n", max_complexity);

	// = NR_REGIONS [ 8 + 4[ #READS + 2[#READS*2READ_LEN] + 2#HAP ]   + 1[ #READS * READ_LEN   +   #HAP * HAP_LEN ] ]
	//int bytes = current_region * (8 + max_reads_nb*(4 + max_read_len) + max_hap_nb*(8 + max_hap_len));
	//printf("Memory to allocate:  %d * (8 + %d*(4 + %d) + %d*(8 + %d))\n", current_region, max_reads_nb , max_read_len, max_hap_nb , max_hap_len);
	
	fclose(stat_file);
	return 1;
}



int add_haplotype(FILE* file) {
	char buffer[BUFFER_SIZE];
	assert(fgets(buffer, BUFFER_SIZE, file));
	char* hap_str = strtok(buffer, ",");
	int hap_length = strlen(hap_str);
	return hap_length;
}




int add_read(FILE* file) {
	char buffer[BUFFER_SIZE];
	assert(fgets(buffer, BUFFER_SIZE, file));
	char* token = strtok(buffer, ",");
	int read_length = strlen(token);
	assert(fgets(buffer, BUFFER_SIZE, file));
	return read_length;
}

