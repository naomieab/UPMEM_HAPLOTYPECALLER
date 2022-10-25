#include "parser.h"
#include "constants.h"


uint32_t offset[NR_REGIONS + 1][OFFSET_SIZE];
uint32_t dpu_region_start_index[NUMBER_DPUS];


uint64_t nr_haplotypes[NUMBER_DPUS]; //an array keeping number of haplotypes in all regions
uint64_t nr_reads[NUMBER_DPUS]; //idem as haplotypes


//TOTAL_READS= maximum total number of reads in a chunk of NUMBER_DPUS regions + MAX_NUMBER OF READS IN ONE REGION  
uint64_t reads_len[TOTAL_READS];
//TOTAL_READS_SIZE = TOTAL_READS * MAX READ LEN  
char reads_array[TOTAL_READS * MAX_READ_LENGTH];
uint64_t haplotypes_len[TOTAL_HAPS];
int64_t haplotypes_val[TOTAL_HAPS];
char haplotypes_array[TOTAL_HAPS * MAX_HAPLOTYPE_LENGTH];
uint32_t priors[TOTAL_READS * MAX_READ_LENGTH * 2];

uint32_t haplotype_region_starts[NUMBER_DPUS][MAX_REGIONS_PER_DPU + 1];
uint32_t read_region_starts[NUMBER_DPUS][MAX_REGIONS_PER_DPU + 1];

int32_t matchToIndel[TOTAL_READS * MAX_READ_LENGTH];



//Addition to keep last region
int tmp_nr_reads;
int tmp_nr_haplotypes;
int tmp_reads_len[MAX_READ_NUM];
char tmp_reads_array[MAX_READ_NUM * MAX_READ_LENGTH];
uint64_t tmp_haplotypes_len[MAX_HAPLOTYPE_NUM];
int64_t tmp_haplotypes_val[MAX_HAPLOTYPE_NUM];
int32_t tmp_priors[2 * MAX_READ_NUM * MAX_READ_LENGTH];
int32_t tmp_matchToIndel[MAX_READ_NUM * MAX_READ_LENGTH];
char tmp_haplotypes_array[MAX_HAPLOTYPE_NUM * MAX_HAPLOTYPE_LENGTH];





//hap_idx is the partial sum of the all the reads in the regions before
//hap_idx+index is the offset where to write this current haplotype
void add_haplotype(FILE* file, int hap_idx, int index) {
	if (hap_idx + index > TOTAL_HAPS) { fprintf(stderr, "Error, number of haps %d is bigger than allocated with TOTAL_HAPS\n", hap_idx + index); }
	char buffer[BUFFER_SIZE];
	assert(fgets(buffer, BUFFER_SIZE, file));

	char* hap_str = strtok(buffer, ",");
	int hap_length = (strlen(hap_str) > MAX_HAPLOTYPE_LENGTH) ? MAX_HAPLOTYPE_LENGTH : strlen(hap_str); //strlen(hap_str);
	haplotypes_val[hap_idx + index] = (int)(log10(1.0 / (hap_length - 1)) * ONE);
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

	assert(fgets(buffer, BUFFER_SIZE, file));
	token = strtok(buffer, ",");
	j = 0;

	for (; token != NULL && j < read_length; token = strtok(NULL, ","), j++) {
		int quality = atoi(token);
		matchToIndel[(read_idx + index) * MAX_READ_LENGTH + j] = (int)(log10(pow((double)10, -(double)quality / 10.0)) * ONE);
	}
}




FILE* read_data(FILE* file, int nr_dpus, int* processed_regions) {
	if (!file) {
		printf("Wrong input file\n");
		return NULL;
	}

	int current_region = 0;
	int current_dpu = 0;
	int hap_idx = 0, read_idx = 0;
	char buffer[BUFFER_SIZE];
	int current_dpu_total_regions = 0;
	int current_dpu_total_complexity = 0;
	dpu_region_start_index[0] = 0;
	nr_reads[0] = 0;
	nr_haplotypes[0] = 0;
	haplotype_region_starts[0][0] = 0;
	read_region_starts[0][0] = 0;
	offset[0][READS_LEN_ARRAY] = 0;
	offset[0][HAPLOTYPES_LEN_VAL_ARRAY] = 0;
	offset[0][READS_ARR] = 0;
	offset[0][HAPS_ARR] = 0;
	offset[0][PRIOR_ARR] = 0;
	int nr_reads_current_region;
	int nr_haplotypes_current_region;
	int region_complexity;
	long int complexity_miss_sources_val[5] = {0,0,0,0,0};
	long int complexity_miss_sources_num[5] = {0,0,0,0,0};
    long int total_complexity = 0;
	long int max_dpu_complexity = 0;
	int sum_read_lengths = 0;
	int sum_haplotypes_lengths = 0;
	printf("There is %d dpus\n", nr_dpus);

	if (*processed_regions != 0 && tmp_nr_haplotypes > 0) {//not first iteration
		nr_haplotypes[0] = tmp_nr_haplotypes;

		nr_reads[0] = tmp_nr_reads;
		sum_haplotypes_lengths = 0;
		sum_read_lengths = 0;
		hap_idx = tmp_nr_haplotypes;
		read_idx = tmp_nr_reads;
		printf("PB1\n");
		for (int i = 0; i < nr_haplotypes[0]; i++) {
			haplotypes_val[i] = tmp_haplotypes_val[i];
			haplotypes_len[i] = tmp_haplotypes_len[i];
			strncpy(&haplotypes_array[i * MAX_HAPLOTYPE_LENGTH], &tmp_haplotypes_array[i * MAX_HAPLOTYPE_LENGTH], MAX_HAPLOTYPE_LENGTH);
			sum_haplotypes_lengths += tmp_haplotypes_len[i];
			printf("PB2\n");
		}
		for (int i = 0; i < nr_reads[0]; i++) {
			reads_len[i] = tmp_reads_len[i];
			strncpy(&reads_array[i * MAX_READ_LENGTH], &tmp_reads_array[i * MAX_READ_LENGTH], MAX_READ_LENGTH);

			for (int j = 0; j < tmp_reads_len[i]; j++) {
				printf("PB4\n");
				priors[i * 2 * MAX_READ_LENGTH + 2 * j] = tmp_priors[i * 2 * MAX_READ_LENGTH + 2 * j];
				priors[i * 2 * MAX_READ_LENGTH + 2 * j + 1] = tmp_priors[i * 2 * MAX_READ_LENGTH + 2 * j + 1];
				matchToIndel[i * MAX_READ_LENGTH + j] = tmp_matchToIndel[i * MAX_READ_LENGTH + j];
			}

			sum_read_lengths += reads_len[i];
		}
		//printf("PB4\n");
		nr_haplotypes[current_dpu] = hap_idx;
		nr_reads[current_dpu] = read_idx;
		current_dpu_total_regions = 1;
		haplotype_region_starts[0][0] = 0;
		haplotype_region_starts[0][1] = nr_haplotypes[0];
		read_region_starts[0][0] = 0;
		read_region_starts[0][1] = nr_reads[0];
		current_dpu_total_complexity = sum_haplotypes_lengths * sum_read_lengths;
		current_region = 1;
		offset[current_region][READS_LEN_ARRAY] = offset[current_region - 1][READS_LEN_ARRAY] + tmp_nr_reads;
		offset[current_region][HAPLOTYPES_LEN_VAL_ARRAY] = offset[current_region - 1][HAPLOTYPES_LEN_VAL_ARRAY] + tmp_nr_haplotypes; //+ (nr_haplotypes_current_region%2==1);
		offset[current_region][READS_ARR] = offset[current_region - 1][READS_ARR] + (tmp_nr_reads * MAX_READ_LENGTH);
		offset[current_region][HAPS_ARR] = offset[current_region - 1][HAPS_ARR] + (tmp_nr_haplotypes * MAX_HAPLOTYPE_LENGTH);
		offset[current_region][PRIOR_ARR] = offset[current_region - 1][PRIOR_ARR] + (tmp_nr_reads * MAX_READ_LENGTH * 2);
	}

	while (current_dpu < nr_dpus && current_region < NR_REGIONS && fgets(buffer, BUFFER_SIZE, file) != 0) {
		/*if (*nr_regions != 0) {
			current_region = 1;

		}*/
		printf("current dpu %d\n", current_dpu);

		sum_read_lengths = 0;
		sum_haplotypes_lengths = 0;
		nr_haplotypes_current_region = atoi(buffer);
		//printf("Number of haps = %s", buffer);
		for (int i = 0; i < nr_haplotypes_current_region; i++) {
			add_haplotype(file, hap_idx, i);
			sum_haplotypes_lengths += haplotypes_len[hap_idx + i];
		}
		hap_idx += nr_haplotypes_current_region;
		assert(nr_haplotypes_current_region <= NR_WRAM_HAPLOTYPES);
		nr_reads_current_region = atoi(fgets(buffer, BUFFER_SIZE, file));

		for (int i = 0; i < nr_reads_current_region; i++) {
			add_read(file, read_idx, i);
			sum_read_lengths += reads_len[read_idx + i];
		}
		read_idx += nr_reads_current_region;
		region_complexity = sum_haplotypes_lengths * sum_read_lengths;
		if (nr_haplotypes_current_region + nr_haplotypes[current_dpu] < MAX_HAPLOTYPE_NUM &&
			nr_reads_current_region + nr_reads[current_dpu] < MAX_READ_NUM &&
			current_dpu_total_regions < MAX_REGIONS_PER_DPU &&
			(current_dpu_total_complexity + region_complexity < TARGET_COMPLEXITY ||
				current_dpu_total_regions == 0)
			) {
			current_dpu_total_complexity += region_complexity;
			nr_haplotypes[current_dpu] += nr_haplotypes_current_region;
			nr_reads[current_dpu] += nr_reads_current_region;
			current_dpu_total_regions++;
			haplotype_region_starts[current_dpu][current_dpu_total_regions] = nr_haplotypes[current_dpu];
			read_region_starts[current_dpu][current_dpu_total_regions] = nr_reads[current_dpu];
		} else {
            total_complexity += current_dpu_total_complexity;
			if (current_dpu_total_complexity+region_complexity >= TARGET_COMPLEXITY) {
				complexity_miss_sources_val[0] += current_dpu_total_complexity;
                complexity_miss_sources_num[0]++;
				printf("\033[32m");
			}
			printf("complexity: %d\t\033[9m%d\033[0m\n", current_dpu_total_complexity, current_dpu_total_complexity + region_complexity);
			if (nr_haplotypes[current_dpu] + nr_haplotypes_current_region >= MAX_HAPLOTYPE_NUM) {
				complexity_miss_sources_val[1] += current_dpu_total_complexity;
                complexity_miss_sources_num[1]++;
				printf("\033[32m");
			}
			printf("nr hapls: %d\t\033[9m%d\033[0m\n", nr_haplotypes[current_dpu], nr_haplotypes[current_dpu] + nr_haplotypes_current_region);
			if (nr_reads[current_dpu] + nr_reads_current_region >= MAX_READ_NUM) {
				complexity_miss_sources_val[2] += current_dpu_total_complexity;
                complexity_miss_sources_num[2]++;
				printf("\033[32m");
			}
			printf("nr reads: %d\t\033[9m%d\033[0m\n", nr_reads[current_dpu], nr_reads[current_dpu] + nr_reads_current_region);
			if (current_dpu_total_regions >= MAX_REGIONS_PER_DPU) {
				complexity_miss_sources_val[3] += current_dpu_total_complexity;
                complexity_miss_sources_num[3]++;
				printf("\033[32m");
			}
			if (current_dpu_total_complexity > max_dpu_complexity) {
				max_dpu_complexity = current_dpu_total_complexity;
			}
			printf("nr regions: %d\t\033[9m%d\033[0m\n\n", current_dpu_total_regions, current_dpu_total_regions + 1);
			current_dpu++;
			dpu_region_start_index[current_dpu] = dpu_region_start_index[current_dpu - 1] + current_dpu_total_regions;
			current_dpu_total_regions = 1;
			//current_dpu_total_reads = nr_reads_current_region;
			current_dpu_total_complexity = region_complexity;
			nr_haplotypes[current_dpu] = nr_haplotypes_current_region;
			nr_reads[current_dpu] = nr_reads_current_region;
			haplotype_region_starts[current_dpu][0] = 0;
			haplotype_region_starts[current_dpu][1] = nr_haplotypes_current_region;
			read_region_starts[current_dpu][0] = 0;
			read_region_starts[current_dpu][1] = nr_reads_current_region;
		}
		current_region++;
		offset[current_region][READS_LEN_ARRAY] = offset[current_region - 1][READS_LEN_ARRAY] + nr_reads_current_region;
		offset[current_region][HAPLOTYPES_LEN_VAL_ARRAY] = offset[current_region - 1][HAPLOTYPES_LEN_VAL_ARRAY] + nr_haplotypes_current_region; //+ (nr_haplotypes_current_region%2==1);
		offset[current_region][READS_ARR] = offset[current_region - 1][READS_ARR] + (nr_reads_current_region * MAX_READ_LENGTH);
		offset[current_region][HAPS_ARR] = offset[current_region - 1][HAPS_ARR] + (nr_haplotypes_current_region * MAX_HAPLOTYPE_LENGTH);
		offset[current_region][PRIOR_ARR] = offset[current_region - 1][PRIOR_ARR] + (nr_reads_current_region * MAX_READ_LENGTH * 2);
	}
	if (current_dpu < nr_dpus && current_dpu_total_complexity > 0) {
        total_complexity += current_dpu_total_complexity;
		printf("complexity: %d\n", current_dpu_total_complexity);
		printf("nr hapls: %d\n", nr_haplotypes[current_dpu]);
		printf("nr reads: %d\n", nr_reads[current_dpu]);
		printf("nr regions: %d\n\n", current_dpu_total_regions);
		complexity_miss_sources_val[4] += current_dpu_total_complexity;
		complexity_miss_sources_num[4]++;
		current_dpu++;
		dpu_region_start_index[current_dpu] = dpu_region_start_index[current_dpu - 1] + current_dpu_total_regions;
		current_dpu_total_regions = 1;
		//current_dpu_total_reads = nr_reads_current_region;
		current_dpu_total_complexity = region_complexity;
		nr_haplotypes[current_dpu] = nr_haplotypes_current_region;
		nr_reads[current_dpu] = nr_reads_current_region;
		haplotype_region_starts[current_dpu][0] = 0;
		haplotype_region_starts[current_dpu][1] = nr_haplotypes_current_region;
		read_region_starts[current_dpu][0] = 0;
		read_region_starts[current_dpu][1] = nr_reads_current_region;
	}

    if (current_dpu>0) {
	    printf("expected wasted dpu time : ~%f%\n", 100*(float)(max_dpu_complexity*current_dpu - total_complexity) / (float)((long int)max_dpu_complexity*current_dpu));
	    printf("loss of evenness due to : \n");
	    printf("complexity : %f%\n", 100*(float)(complexity_miss_sources_num[0]*max_dpu_complexity - complexity_miss_sources_val[0])/(float)((long int)max_dpu_complexity*current_dpu));
	    printf("nr hapls   : %f%\n", 100*(float)(complexity_miss_sources_num[1]*max_dpu_complexity - complexity_miss_sources_val[1])/(float)((long int)max_dpu_complexity*current_dpu));
	    printf("nr reads   : %f%\n", 100*(float)(complexity_miss_sources_num[2]*max_dpu_complexity - complexity_miss_sources_val[2])/(float)((long int)max_dpu_complexity*current_dpu));
	    printf("nr regions : %f%\n", 100*(float)(complexity_miss_sources_num[3]*max_dpu_complexity - complexity_miss_sources_val[3])/(float)((long int)max_dpu_complexity*current_dpu));
	    printf("last dpu   : %f%\n", 100*(float)(complexity_miss_sources_num[4]*max_dpu_complexity - complexity_miss_sources_val[4])/(float)((long int)max_dpu_complexity*current_dpu));
	    printf("highest dpu complexity : %d\n", max_dpu_complexity);
	    printf("Dpus used  : %d/%d\n", current_dpu, nr_dpus);
        printf("regions sent : %d\n", current_region);
    }

	//cover extreme case of last region
    if (current_dpu == nr_dpus) {// Only if all dpus were used
        tmp_nr_haplotypes = nr_haplotypes_current_region;
        tmp_nr_reads = nr_reads_current_region;
        for (int i = 0; i < tmp_nr_haplotypes; i++) {
            tmp_haplotypes_val[i] = haplotypes_val[hap_idx - tmp_nr_haplotypes + i];
            tmp_haplotypes_len[i] = haplotypes_len[hap_idx - tmp_nr_haplotypes + i];
            strncpy(&tmp_haplotypes_array[i * MAX_HAPLOTYPE_LENGTH], &haplotypes_array[(hap_idx - tmp_nr_haplotypes + i) * MAX_HAPLOTYPE_LENGTH], MAX_HAPLOTYPE_LENGTH);
        }
        for (int i = 0; i < tmp_nr_reads; i++) {
            tmp_reads_len[i] = reads_len[read_idx - tmp_nr_reads + i];
            strncpy(&tmp_reads_array[i * MAX_READ_LENGTH], &reads_array[(read_idx - tmp_nr_reads + i) * MAX_READ_LENGTH], MAX_READ_LENGTH);
            for (int j = 0; j < tmp_reads_len[i]; j++) {
                tmp_priors[i * 2 * MAX_READ_LENGTH + 2 * j] = priors[(read_idx - tmp_nr_reads + i) * 2 * MAX_READ_LENGTH + 2 * j];
                tmp_priors[i * 2 * MAX_READ_LENGTH + 2 * j + 1] = priors[(read_idx - tmp_nr_reads + i) * 2 * MAX_READ_LENGTH + 2 * j + 1];
                tmp_matchToIndel[i * MAX_READ_LENGTH + j] = matchToIndel[(read_idx - tmp_nr_reads + i) * MAX_READ_LENGTH + j];
            }
        }
        // printf("last region saved for next dpu run\n");
    } else {
        tmp_nr_haplotypes = 0;
    }


	// Set to 0 the sizes of all unused regions
	while (current_dpu < NUMBER_DPUS) {
		nr_reads[current_dpu] = 0;
		nr_haplotypes[current_dpu] = 0;
		current_dpu++;
	}
	*processed_regions = *processed_regions + current_region;
	printf("Sent %d regions, max complexity = %d\n", current_region, max_dpu_complexity);
	return file;
}




void free_mem(FILE* file) {
	fclose(file);
}
