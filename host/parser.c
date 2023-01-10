#include "parser.h"
#include "constants.h"
#include "log.h"
#include "buffers.h"
#include <limits.h>
#include <float.h>
#define HAPLOTYPE_STEP 3
#define min(X,Y) (((X) < (Y)) ? (X) : (Y))

//data of scanned region
uint32_t nr_reads_region;
uint32_t nr_haplotypes_region;
char reads_lines[2 * MAX_READS_REGION][BUFFER_SIZE];
uint32_t reads_len[MAX_READS_REGION];
//uint32_t haps_len[MAX_HAPLOTYPE_REGION];
char haplotypes_lines[MAX_HAPLOTYPE_REGION][BUFFER_SIZE];
int priors_arr[MAX_READS_REGION][MAX_READ_LENGTH];
int transitions_arr[MAX_READS_REGION][MAX_READ_LENGTH];
uint32_t line;
/*
* Read a region from input file
* @return complexity of the region
*/
FILE* scan_region(FILE* file, long long int* region_complexity, int current_region) {
	int total_read_lengths = 0, total_hap_lengths = 0;
	char buffer[BUFFER_SIZE];
	fgets(buffer, BUFFER_SIZE, file);
	line++;
	nr_haplotypes_region = atoi(buffer);
	for (int j = 0; j < nr_haplotypes_region; j++) { //scan haplotypes
		fgets(buffer, BUFFER_SIZE, file);
		line++;
		strcpy(haplotypes_lines[j], buffer);
		total_hap_lengths += strlen(buffer) - 1;
	}
	fgets(buffer, BUFFER_SIZE, file);
	line++;
	nr_reads_region = atoi(buffer);
	assert(nr_reads_region <= MAX_READS_REGION);
	for (int j = 0; j < nr_reads_region; j++) { //scan reads
		fgets(buffer, BUFFER_SIZE, file);
		line++;
		strcpy(reads_lines[2 * j], buffer);
		char* tmp = strchr(buffer, ',');
		reads_len[j] = tmp - buffer;//tmp is the adress in buffer of the first occurrence of a comma => the difference here is the length of the read
		total_read_lengths += reads_len[j];


		//extract qualities
		char* comma = strchr(tmp + 1, ',');
		char* previous = tmp + 1;
		int index_in_str = comma - buffer + 2;
		int length_of_substring = comma - previous;
		if (length_of_substring >= 4) {
			LOG_INFO(">4: region %d\n", current_region);
		}
		char str_val[4];
		int k = 0, quality;
		assert(length_of_substring < 4);
		while (comma != NULL) {
			strncpy(str_val, previous, length_of_substring);
			str_val[length_of_substring] = '\0';
			priors_arr[j][k] = atoi(str_val);
			previous = comma + 2;//+1 because comma is the current comma, and +1 because there is a space between each
			comma = strchr(comma + 2, ',');
			length_of_substring = comma - previous;
			k++;
		}
		//last one (doesn't have comma after that)
		str_val[0] = *previous;
		str_val[1] = '\0';
		//TODO: check if not need to change the condition to \n
		if (*(previous + 1) != '\0') { str_val[1] = *(previous + 1); str_val[2] = '\0'; }
		priors_arr[j][k] = atoi(str_val);

		fgets(buffer, BUFFER_SIZE, file);
		line++;
		strcpy(reads_lines[2 * j + 1], buffer);

		//extract transitions
		previous = buffer;
		comma = buffer;
		k = 0;

		//TODO: reorganize inside the loop to be able to remove the if condition 
		while (comma != NULL) {
			comma = strchr(comma + 2, ',');
			if (comma != NULL) {
				length_of_substring = comma - previous;
				assert(length_of_substring < 3);
				strncpy(str_val, previous, length_of_substring);
				str_val[length_of_substring] = '\0';
				previous = comma + 2;
				transitions_arr[j][k] = atoi(str_val);
				k++;
			}
		}
		transitions_arr[j][k] = 45;
	}
	*region_complexity = (long long int) total_read_lengths * (long long int) total_hap_lengths;
	return file;
}

/*
* This function loads the relative informations to send a part/whole region
* read/hap_idx are the index of the read/hap in the buffers
* read/hap_nb are the number of reads/haps to send
* read/hap_offset is the offset of reads/hap in the array to load
*/
void send_region(int current_dpu, int read_idx, int hap_idx, int read_nb, int hap_nb, int read_offset, int hap_offset, int region_idx, int region_sub_idx, bool last_subregion) {
	assert(hap_nb < NR_WRAM_HAPLOTYPES);
	int dpu_region = dpu_regions_buffer[current_dpu].nr_regions;
	dpu_regions_buffer[current_dpu].region_shapes[dpu_region].region_index = region_idx;
	dpu_regions_buffer[current_dpu].region_shapes[dpu_region].nr_reads = read_nb;
	dpu_regions_buffer[current_dpu].region_shapes[dpu_region].nr_haplotypes = hap_nb;
	dpu_regions_buffer[current_dpu].region_shapes[dpu_region].read_offset = read_offset;
	dpu_regions_buffer[current_dpu].region_shapes[dpu_region].hapl_offset = hap_offset;
	dpu_regions_buffer[current_dpu].region_shapes[dpu_region].total_nr_subregions = region_sub_idx; //the sub region index for current region
	dpu_regions_buffer[current_dpu].region_shapes[dpu_region].last_subregion = last_subregion;
	dpu_regions_buffer[current_dpu].region_shapes[dpu_region].total_reads_region = nr_reads_region;
	dpu_regions_buffer[current_dpu].region_shapes[dpu_region].total_haps_region = nr_haplotypes_region;

	dpu_regions_buffer[current_dpu].haplotype_region_starts[dpu_region] = dpu_regions_buffer[current_dpu].nr_haplotypes - hap_nb; //substract because we have already updated the nb of hap in dpu
	dpu_regions_buffer[current_dpu].read_region_starts[dpu_region] = dpu_regions_buffer[current_dpu].nr_reads - read_nb;

	if (dpu_region < MAX_REGIONS_PER_DPU) {//TODO:the condition can be removed
		dpu_regions_buffer[current_dpu].haplotype_region_starts[dpu_region + 1] = dpu_regions_buffer[current_dpu].nr_haplotypes;
		dpu_regions_buffer[current_dpu].read_region_starts[dpu_region + 1] = dpu_regions_buffer[current_dpu].nr_reads;
	}
	assert(dpu_regions_buffer[current_dpu].haplotype_region_starts[0] == 0);
	assert(dpu_regions_buffer[current_dpu].read_region_starts[0] == 0);
	//fill haplotypes buffers
	for (int i = 0; i < hap_nb; i++) {
		int hap_length = strlen(haplotypes_lines[hap_offset + i]);//strlen(strtok(haplotypes_lines[hap_offset + i], ","));
		haplotypes_len_buffer[hap_idx + i] = hap_length;
		assert(hap_length <= MAX_HAPLOTYPE_LENGTH);
		// LOG_DEBUG("len:%d %s\n", hap_length, haplotypes_lines[hap_offset + i]);
		haplotypes_val_buffer[hap_idx + i] = (int64_t)(log10(1.0 / (hap_length - 1)) * ONE);
		strncpy(&haplotypes_array_buffer[(hap_idx + i) * MAX_HAPLOTYPE_LENGTH], haplotypes_lines[hap_offset + i], MAX_HAPLOTYPE_LENGTH);
	}

	if (dpu_regions_buffer[current_dpu].nr_regions == 0) {
		assert(hap_idx + MAX_HAPLOTYPE_NUM < TOTAL_HAPS);
		dpu_regions_buffer[current_dpu].haplotypes_len = &haplotypes_len_buffer[hap_idx];
		dpu_regions_buffer[current_dpu].haplotypes_val = &haplotypes_val_buffer[hap_idx];
		dpu_regions_buffer[current_dpu].haplotypes_array = &haplotypes_array_buffer[hap_idx * MAX_HAPLOTYPE_LENGTH];
	}

	//fill reads buffers
	for (int i = 0; i < read_nb; i++) {
		int read_length = reads_len[read_offset + i];
		reads_len_buffer[read_idx + i] = read_length;

		strncpy(&reads_array_buffer[(read_idx + i) * MAX_READ_LENGTH], reads_lines[2 * (read_offset + i)], read_length);
		reads_array_buffer[(read_idx + i) * MAX_READ_LENGTH + read_length] = '\0';

		// LOG_DEBUG("Read %d:%s", i, reads_lines[2 * (read_offset + i)]);

		double probLog10, errorProbLog10;
		int quality;
		for (int j = 0; j < read_length; j++) {
			quality = priors_arr[read_offset + i][j];
			probLog10 = log10(1 - pow((double)10, -(double)quality / 10.0));
			errorProbLog10 = log10(pow((double)10, -(double)quality / 10.0)) - log10(3);
			priors[(read_idx + i) * 2 * MAX_READ_LENGTH + 2 * j] = (int)(probLog10 * ONE);
			priors[(read_idx + i) * 2 * MAX_READ_LENGTH + 2 * j + 1] = (int)(errorProbLog10 * ONE);
			match_to_indel_buffer[(read_idx + i) * MAX_READ_LENGTH + j] = (int)(log10(pow((double)10, -(double)transitions_arr[read_offset + i][j] / 10.0)) * ONE);
		}
	}
	if (dpu_regions_buffer[current_dpu].nr_regions == 0) {
		dpu_regions_buffer[current_dpu].reads_len = &reads_len_buffer[read_idx];
		dpu_regions_buffer[current_dpu].reads_array = &reads_array_buffer[read_idx * MAX_READ_LENGTH];
		dpu_regions_buffer[current_dpu].priors = &priors[read_idx * 2 * MAX_READ_LENGTH];
		dpu_regions_buffer[current_dpu].match_to_indel = &match_to_indel_buffer[read_idx * MAX_READ_LENGTH];
	}

	return;
}

void* read_data(void* input_file) {
	FILE* file = (FILE*)input_file;
	int current_region = -1;
	long long int region_complexity;
	int current_dpu_left_complexity = TARGET_COMPLEXITY;
	int read_idx = 0, hap_idx = 0;
	int read_offset = 0, hap_offset = 0;
	line = 0;
	//initialize queue
	queue_init(&dpu_regions_queue, DPU_INPUT_BUFFER_SIZE);

	int current_dpu = queue_put(&dpu_regions_queue);

	//initialization for first dpu
	dpu_regions_buffer[current_dpu].dpu_inactive = 0;
	dpu_regions_buffer[current_dpu].first_region_index = 0;
	dpu_regions_buffer[current_dpu].nr_regions = 0;
	dpu_regions_buffer[current_dpu].nr_reads = 0;
	dpu_regions_buffer[current_dpu].nr_haplotypes = 0;

	while (current_region < TOTAL_REGIONS - 1) {//-1 because we do the incrementation at the begining of the loop
		file = scan_region(file, &region_complexity, current_region);
		hap_offset = 0;
		read_offset = 0;
		current_region++;

		if (region_complexity <= current_dpu_left_complexity &&
			dpu_regions_buffer[current_dpu].nr_reads + nr_reads_region <= MAX_READ_NUM &&
			dpu_regions_buffer[current_dpu].nr_haplotypes + nr_haplotypes_region <= MAX_HAPLOTYPE_NUM &&
			dpu_regions_buffer[current_dpu].nr_regions + 1 <= MAX_REGIONS_PER_DPU &&
			nr_haplotypes_region < NR_WRAM_HAPLOTYPES) {

			current_dpu_left_complexity -= region_complexity;
			dpu_regions_buffer[current_dpu].nr_reads += nr_reads_region;
			dpu_regions_buffer[current_dpu].nr_haplotypes += nr_haplotypes_region;
			send_region(current_dpu, read_idx, hap_idx, nr_reads_region, nr_haplotypes_region, 0, 0, current_region, 1, 1);

			dpu_regions_buffer[current_dpu].nr_regions++;

			read_idx += nr_reads_region;
			hap_idx += nr_haplotypes_region;
		}
		else { //partage et si ca rentre plus chsange de dpu
		   //tant que tu peux en rentrer rajoute
		   //pour garder la meme function on fait dabord le passage ici pour decider cb de read et hap ( a checker niveau temps si cest trop faut l'inclure dans la fonction)
			int current_sub_region = 1;
			int tmp_hap = 0;
			while (region_complexity > 0) {
				int hap_to_add = min(nr_haplotypes_region - hap_offset, HAPLOTYPE_STEP); //first item is the number of remaining haps in region
				if (dpu_regions_buffer[current_dpu].nr_haplotypes + ((tmp_hap != 0) ? tmp_hap : hap_to_add) > MAX_HAPLOTYPE_NUM ||
					dpu_regions_buffer[current_dpu].nr_reads >= MAX_READ_NUM ||
					dpu_regions_buffer[current_dpu].nr_regions + 1 > MAX_REGIONS_PER_DPU ||
					current_dpu_left_complexity < 0) {
					uint64_t* pnt = dpu_regions_buffer[current_dpu].reads_len;
					for (int k = 0; k < dpu_regions_buffer[current_dpu].nr_reads; k++) {
						assert(*pnt != 0);
						pnt++;
					}
					pnt = dpu_regions_buffer[current_dpu].haplotypes_len;
					for (int k = 0; k < dpu_regions_buffer[current_dpu].nr_haplotypes; k++) {
						assert(*pnt != 0);
						pnt += 1;
					}
					queue_make_available(&dpu_regions_queue, current_dpu);
					LOG_DEBUG("NewDPU\n");
					current_dpu = queue_put(&dpu_regions_queue);
					current_dpu_left_complexity = TARGET_COMPLEXITY;
					dpu_regions_buffer[current_dpu].nr_reads = 0;
					dpu_regions_buffer[current_dpu].nr_haplotypes = 0;
					dpu_regions_buffer[current_dpu].nr_regions = 0;
					dpu_regions_buffer[current_dpu].first_region_index = current_region;
					dpu_regions_buffer[current_dpu].dpu_inactive = 0;
					if (read_idx + nr_reads_region > TOTAL_READS - MAX_READ_NUM || hap_idx + nr_haplotypes_region > TOTAL_HAPS - MAX_HAPLOTYPE_NUM) {//
						read_idx = 0;
						hap_idx = 0;
					}
				}
				int partial_read_sum = 0;
				int read_cnt = 0, hap_cnt = 0;
				while (current_dpu_left_complexity > 0 &&
					dpu_regions_buffer[current_dpu].nr_reads < MAX_READ_NUM &&
					(hap_offset < nr_haplotypes_region || read_offset < nr_reads_region) &&
					dpu_regions_buffer[current_dpu].nr_haplotypes + ((tmp_hap != 0) ? tmp_hap : hap_to_add) <= MAX_HAPLOTYPE_NUM &&
					dpu_regions_buffer[current_dpu].nr_regions + 1 <= MAX_REGIONS_PER_DPU) {
					if (tmp_hap != 0) { hap_to_add = tmp_hap; }
					else {
						hap_to_add = min(nr_haplotypes_region - hap_offset, HAPLOTYPE_STEP);
					}
					partial_read_sum = 0;
					dpu_regions_buffer[current_dpu].nr_haplotypes += hap_to_add;

					int length_hap = 0;
					for (int i = 0; i < hap_to_add; i++) {
						length_hap += strlen(haplotypes_lines[hap_offset + i]) - 1;
					}
					hap_cnt += hap_to_add;
					int reads_length_goal = current_dpu_left_complexity / length_hap;

					while (partial_read_sum <= reads_length_goal && read_offset + read_cnt < nr_reads_region &&
						dpu_regions_buffer[current_dpu].nr_reads + read_cnt < MAX_READ_NUM) { //add reads

						partial_read_sum += reads_len[read_offset + read_cnt];
						read_cnt++;
						dpu_regions_buffer[current_dpu].nr_reads++;
					}

					while (length_hap * partial_read_sum < current_dpu_left_complexity &&
						hap_offset + hap_cnt < nr_haplotypes_region &&
						dpu_regions_buffer[current_dpu].nr_haplotypes + hap_cnt < MAX_HAPLOTYPE_NUM &&
						read_offset == 0 &&
						hap_cnt < NR_WRAM_HAPLOTYPES - 1) {//check if we can add haplotypes

						length_hap += strlen(haplotypes_lines[hap_offset + hap_cnt]) - 1;
						hap_cnt++;
						dpu_regions_buffer[current_dpu].nr_haplotypes++;

					}
					current_dpu_left_complexity -= (length_hap * partial_read_sum);
					region_complexity -= (length_hap * partial_read_sum);
					send_region(current_dpu, read_idx, hap_idx, read_cnt, hap_cnt, read_offset, hap_offset, current_region, current_sub_region++, (region_complexity <= 0));
					//DEBUG PURPOSE
          /*if (current_region == 7173 || current_region == 7174) {
						LOG_INFO("Send region on idx %d read %d hap with %d read, %d hap, %d read_offset, %d hap_offset,  region_idx=%d, int region_sub_idx=%d\n", read_idx, hap_idx, read_cnt, hap_cnt, read_offset, hap_offset, current_region, current_sub_region);
						for (int p = 0; p < hap_cnt; p++) {
							LOG_INFO("%d, %s", strlen(dpu_regions_buffer[current_dpu].haplotypes_array + p * MAX_HAPLOTYPE_LENGTH), (dpu_regions_buffer[current_dpu].haplotypes_array + p * MAX_HAPLOTYPE_LENGTH));
						}
					}*/
					dpu_regions_buffer[current_dpu].nr_regions++;

					read_idx += read_cnt;
					hap_idx += hap_cnt;


					if (hap_cnt != HAPLOTYPE_STEP && read_offset + read_cnt != nr_reads_region) {
						tmp_hap = hap_cnt;
					}
					else { tmp_hap = 0; }

					if (read_offset + read_cnt != nr_reads_region) {
						read_offset += read_cnt;
					}
					else if (hap_offset + hap_cnt == nr_haplotypes_region) {
						read_offset += read_cnt;
						hap_offset += hap_cnt;
						read_cnt = 0;
						hap_cnt = 0;
						continue;
					}
					else {
						hap_offset += hap_cnt;
						read_cnt = 0;
						read_offset = 0;
					}
					read_cnt = 0;
					hap_cnt = 0;

				}

			}
		}

	}
	LOG_INFO("\033[42mFinished regions scan\033[0m\n");
    queue_close(&dpu_regions_queue, MAX_RANKS+1);
	return NULL;
}
