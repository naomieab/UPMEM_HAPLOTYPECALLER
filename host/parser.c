#include "parser.h"
#include "constants.h"
#include "log.h"
#include "buffers.h"

#define HAPLOTYPE_STEP 3

//data of scanned region
uint32_t nr_reads_region;
uint32_t nr_haplotypes_region;
char reads_lines[2 * MAX_READS_REGION][BUFFER_SIZE];
char haplotypes_lines[MAX_HAPLOTYPE_REGION][BUFFER_SIZE];

/*
* Read a region from input file 
* @return complexity of the region
*/
FILE* scan_region(FILE* file, int* region_complexity){
	printf("Scan new region\n");
	int total_read_lengths = 0, total_hap_lengths = 0;
	char buffer[BUFFER_SIZE];
	fgets(buffer, BUFFER_SIZE, file);
	nr_haplotypes_region = atoi(buffer);
	printf("Haplotypes:%d  -  ", nr_haplotypes_region);
	assert(nr_haplotypes_region <= MAX_HAPLOTYPE_REGION); 
	for (int j = 0; j < nr_haplotypes_region; j++) { //scan haplotypes
		fgets(buffer, BUFFER_SIZE, file);
		strcpy(haplotypes_lines[j], buffer);
		total_hap_lengths += strlen(strtok(buffer, ","));
		//printf("Total hap length = %d\n", total_hap_lengths);
	}
	fgets(buffer, BUFFER_SIZE, file);
	nr_reads_region = atoi(buffer);
	printf("Reads:%d\n", nr_reads_region);
	assert(nr_reads_region <= MAX_READS_REGION);
	for (int j = 0; j < nr_reads_region; j++) { //scan reads
		fgets(buffer, BUFFER_SIZE, file);
		strcpy(reads_lines[2 * j], buffer);
		total_read_lengths += strlen(strtok(buffer, ","));
		fgets(buffer, BUFFER_SIZE, file);
		strcpy(reads_lines[2 * j + 1], buffer);
	}
	//printf("Total read length %d and total hap length %d\n", total_read_lengths, total_hap_lengths);
	*region_complexity = total_read_lengths * total_hap_lengths;
	return file;
}

/*
* This function loads the relative informations to send a part/whole region
* read/hap_idx are the index of the read/hap in the buffers
* read/hap_nb are the number of reads/haps to send
* read/hap_offset is the offset of reads/hap in the array to load
*/
void send_region(int current_dpu, int read_idx, int hap_idx, int read_nb, int hap_nb, int read_offset, int hap_offset, int region_idx, int region_sub_idx, bool last_subregion) {
	printf("Send region with read_idx=%d, int hap_idx=%d, int read_nb=%d, int hap_nb=%d, int read_offset=%d, int hap_offset=%d, int region_idx=%d, int region_sub_idx=%d\n", read_idx,  hap_idx, read_nb, hap_nb, read_offset, hap_offset, region_idx, region_sub_idx);
	if (last_subregion) { printf("LAST SUBREGION\n\n"); }
	int dpu_region = dpu_regions_buffer[current_dpu].nr_regions;
	struct region_shape_t region_shape = dpu_regions_buffer[current_dpu].region_shapes[dpu_region];

	region_shape.region_index = region_idx;
	region_shape.nr_reads = read_nb;
	region_shape.nr_haplotypes = hap_nb;
	region_shape.read_offset = read_offset;
	region_shape.hapl_offset = hap_offset;
	region_shape.total_nr_subregions = region_sub_idx; //the sub region index for current region
	region_shape.last_subregion = last_subregion;
	region_shape.total_reads_region = nr_reads_region;
	region_shape.total_haps_region = nr_haplotypes_region;

	dpu_regions_buffer[current_dpu].haplotype_region_starts[dpu_region] = dpu_regions_buffer[current_dpu].nr_haplotypes - hap_nb; //substract because we have already updated the nb of hap in dpu
	dpu_regions_buffer[current_dpu].read_region_starts[dpu_region] = dpu_regions_buffer[current_dpu].nr_reads - read_nb;

	//fill haplotypes buffers
	for (int i = 0; i < hap_nb; i++) {
		int hap_length = strlen(strtok(haplotypes_lines[hap_offset + i], ","));
		haplotypes_len_buffer[hap_idx + i] = hap_length;
		haplotypes_val_buffer[hap_idx + i] = (int)(log10(1.0 / (hap_length - 1)) * ONE);
		strncpy(&haplotypes_array_buffer[(hap_idx + i) * MAX_HAPLOTYPE_LENGTH], haplotypes_lines[hap_offset + i], MAX_HAPLOTYPE_LENGTH);
	}

	dpu_regions_buffer[current_dpu].haplotypes_len = &haplotypes_len_buffer[hap_idx];
	dpu_regions_buffer[current_dpu].haplotypes_val = &haplotypes_val_buffer[hap_idx];
	dpu_regions_buffer[current_dpu].haplotypes_array = &haplotypes_array_buffer[hap_idx * MAX_HAPLOTYPE_LENGTH];


	//fill reads buffers
	for (int i = 0; i < read_nb; i++) {
		int read_length = strlen(strtok(reads_lines[2 * (read_offset + i)], ","));
		reads_len_buffer[read_idx + i] = read_length;
		strcpy(&reads_array_buffer[(read_idx + i) * MAX_READ_LENGTH], strtok(reads_lines[2*(read_offset + i)], ","));
		char* prior;
		char* indel;
		int j;
		for (j = 0, prior = strtok(NULL, ","), indel = strtok(reads_lines[2 * (read_offset + i) + 1], ","); prior != NULL && j < read_length; prior = strtok(NULL, ","), indel = strtok(NULL, ","), j++) {
			int quality = atoi(prior);
			int transition = atoi(indel);
			double probLog10 = log10(1 - pow((double)10, -(double)quality / 10.0));
			double errorProbLog10 = log10(pow((double)10, -(double)quality / 10.0)) - log10(3);
			priors[(read_idx + i) * 2 * MAX_READ_LENGTH + 2 * j] = (int)(probLog10 * ONE);
			priors[(read_idx + i) * 2 * MAX_READ_LENGTH + 2 * j + 1] = (int)(errorProbLog10 * ONE);
			match_to_indel_buffer[(read_idx + i) * MAX_READ_LENGTH + j] = (int)(log10(pow((double)10, -(double)transition / 10.0)) * ONE);
		}

	}

	dpu_regions_buffer[current_dpu].reads_len = &reads_len_buffer[read_idx];
	dpu_regions_buffer[current_dpu].reads_array = &reads_array_buffer[read_idx * MAX_READ_LENGTH];
	dpu_regions_buffer[current_dpu].priors = &priors[read_idx * 2 * MAX_READ_LENGTH];
	dpu_regions_buffer[current_dpu].match_to_indel = &match_to_indel_buffer[read_idx * MAX_READ_LENGTH];
	
	return;
}

void read_data(FILE* file, int nr_dpus) {
	int current_region = -1;
	int region_complexity, current_dpu_left_complexity = TARGET_COMPLEXITY;
	int read_idx = 0, hap_idx = 0;
	int read_offset = 0, hap_offset = 0;

	//initialize queue
	queue_init(&dpu_regions_queue, DPU_OUTPUT_BUFFER_SIZE);

	int current_dpu = queue_put(&dpu_regions_queue);

	//initialization for first dpu
	dpu_regions_buffer[current_dpu].dpu_inactive = 0;
	dpu_regions_buffer[current_dpu].first_region_index = 0;
	dpu_regions_buffer[current_dpu].nr_regions = 0;
	dpu_regions_buffer[current_dpu].nr_reads = 0;
	dpu_regions_buffer[current_dpu].nr_haplotypes = 0;

	while (current_region < TOTAL_REGIONS-1) {//-1 because we do the incrementation at the begining of the loop
		file = scan_region(file, &region_complexity);
		hap_offset = 0;
		read_offset = 0;
		current_region++;
		printf("Region complexity is %d\n", region_complexity);
		if (region_complexity <= current_dpu_left_complexity &&
			dpu_regions_buffer[current_dpu].nr_reads + nr_reads_region <= MAX_READ_NUM &&
			dpu_regions_buffer[current_dpu].nr_haplotypes + nr_haplotypes_region <= MAX_HAPLOTYPE_NUM &&
			dpu_regions_buffer[current_dpu].nr_regions + 1 <= MAX_REGIONS_PER_DPU) {

			current_dpu_left_complexity -= region_complexity;
			dpu_regions_buffer[current_dpu].nr_reads += nr_reads_region;
			dpu_regions_buffer[current_dpu].nr_haplotypes += nr_haplotypes_region;

			send_region(current_dpu, read_idx, hap_idx, nr_reads_region, nr_haplotypes_region, 0, 0, current_region, 0, (region_complexity <= 0));
			
			dpu_regions_buffer[current_dpu].nr_regions++;

			read_idx += nr_reads_region;
			hap_idx += nr_haplotypes_region;
		}else{ //partage et si ca rentre plus chsange de dpu
			//tant que tu peux en rentrer rajoute
			//pour garder la meme function on fait dabord le passage ici pour decider cb de read et hap ( a checker niveau temps si cest trop faut l'inclure dans la fonction)
			int current_sub_region = 0;
			while (region_complexity > 0) {
				printf("Enter splitting of region dpu_complexity %d and region complexity %d\n", current_dpu_left_complexity, region_complexity);
				if (dpu_regions_buffer[current_dpu].nr_haplotypes + HAPLOTYPE_STEP > MAX_HAPLOTYPE_NUM ||
					dpu_regions_buffer[current_dpu].nr_reads >= MAX_READ_NUM ||
					dpu_regions_buffer[current_dpu].nr_regions + 1 > MAX_REGIONS_PER_DPU ||
					current_dpu_left_complexity < 0) {
					queue_make_available(&dpu_regions_queue, current_dpu);
					fprintf(stderr, "Allocate a NEW DPU dpu left complexity %d\n", current_dpu_left_complexity);
					current_dpu = queue_put(&dpu_regions_queue);
					//if (current_sub_region != 0) { current_sub_region++; }
					//if (current_dpu == nr_dpus) { current_dpu = 0; } //TODO: must wait for queue to free something?
					current_dpu_left_complexity = TARGET_COMPLEXITY;
					dpu_regions_buffer[current_dpu].nr_reads = 0;
					dpu_regions_buffer[current_dpu].nr_haplotypes = 0;
					dpu_regions_buffer[current_dpu].nr_regions = 0;
					dpu_regions_buffer[current_dpu].first_region_index = current_region;
					dpu_regions_buffer[current_dpu].dpu_inactive = 0;
					if (read_idx + nr_reads_region > TOTAL_READS || hap_idx + nr_haplotypes_region > TOTAL_HAPS) {//
						read_idx = 0;
						hap_idx = 0;
					}
				}
				
				int partial_read_sum = 0;
				int read_cnt = 0, hap_cnt = 0;
				printf("Hap idx %d read idx %d dpu nr reads %d", hap_idx, read_idx, dpu_regions_buffer[current_dpu].nr_reads);
				while (current_dpu_left_complexity > 0 &&
					dpu_regions_buffer[current_dpu].nr_reads < MAX_READ_NUM &&
					(hap_offset < nr_haplotypes_region || read_offset < nr_reads_region)) {
					partial_read_sum = 0;
					dpu_regions_buffer[current_dpu].nr_haplotypes += 3;
					int length_hap = 0;
					for (int i = 0; i < HAPLOTYPE_STEP; i++) {
						length_hap += strlen(strtok(haplotypes_lines[hap_offset + i], ","));
					}
					hap_cnt += HAPLOTYPE_STEP;
					int reads_length_goal = current_dpu_left_complexity / length_hap;
					
					while (partial_read_sum < reads_length_goal && read_offset + read_cnt < nr_reads_region &&
						dpu_regions_buffer[current_dpu].nr_reads + read_cnt < MAX_READ_NUM) { //add reads

						partial_read_sum += strlen(strtok(reads_lines[2 * (read_offset + read_cnt)], ","));
						read_cnt++;
					}
				
					while (length_hap * partial_read_sum < current_dpu_left_complexity &&
						hap_offset + hap_cnt < nr_haplotypes_region &&
						dpu_regions_buffer[current_dpu].nr_haplotypes + hap_cnt < MAX_HAPLOTYPE_NUM &&
						read_offset == 0) {//check if we can add haplotypes

						length_hap += strlen(strtok(haplotypes_lines[hap_offset + hap_cnt], ","));
						hap_cnt++;
					}
					current_dpu_left_complexity -= (length_hap * partial_read_sum);
					region_complexity -= (length_hap * partial_read_sum);

					send_region(current_dpu, read_idx, hap_idx, read_cnt, hap_cnt, read_offset, hap_offset, current_region, current_sub_region++, (region_complexity<=0));
					dpu_regions_buffer[current_dpu].nr_regions++;

					printf("Current dpu left complexity %d\n", current_dpu_left_complexity);
					read_idx += read_cnt;
					hap_idx += hap_cnt;

					
					
					if (read_offset+read_cnt != nr_reads_region) {
						//hap_offset = 0;
						read_offset += read_cnt;
					}
					else if (hap_offset+hap_cnt == nr_haplotypes_region) {//if (read_offset+read_cnt == nr_reads_region && hap offset==nrhaps_region)
						printf("Finish dpu left = %d\n", current_dpu_left_complexity);
						read_offset += read_cnt;
						hap_offset += hap_cnt;
						read_cnt = 0; 
						hap_cnt = 0;
						continue;
					}
					else {// (read_offset == nr_reads_region) {
						hap_offset += hap_cnt;
						read_cnt = 0;
						read_offset = 0;
					}

					printf("Has send region now read cnt = %d and hap cnt = %d\n", read_idx, hap_idx);
					read_cnt = 0;
					hap_cnt = 0;
		
				}

			}
		}
	}
}
