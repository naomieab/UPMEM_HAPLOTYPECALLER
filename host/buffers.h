#include "circular_queue.h"
#include "constants.h"
#include <stdint.h>

struct region_shape_t {
	// The index of this region (or the index of the region wich this subregion is from)
	uint32_t region_index;
	uint32_t nr_reads;
	uint32_t nr_haplotypes;
	uint32_t read_offset;
	uint32_t hapl_offset;
	uint32_t total_nr_subregions;// The number of subregions from the same main region
};

struct dpu_regions_t {
	uint64_t  dpu_inactive;
	uint32_t  first_region_index;
	uint32_t  nr_regions;
	uint64_t  nr_reads;
	uint64_t  nr_haplotypes;
	struct    region_shape_t region_shapes[MAX_REGIONS_PER_DPU];
	uint32_t  haplotype_region_starts[MAX_REGIONS_PER_DPU+1];
	uint32_t  read_region_starts[MAX_REGIONS_PER_DPU+1];
	uint64_t* reads_len;
	char*     reads_array;
	uint64_t* haplotypes_len;
	uint64_t* haplotypes_val;
	char*     haplotypes_array;
	uint32_t* priors;
	int32_t*  match_to_indel;
};

struct dpu_results_t {
	uint32_t first_region_index;
	uint32_t nr_regions;
	struct   region_shape_t region_shapes[MAX_REGIONS_PER_DPU];
	int64_t  likelihoods[MAX_READ_NUM][MAX_HAPLOTYPE_NUM];
	uint64_t nr_cycles;
};

struct dpu_work_t {
	uint32_t first_region_index;
	uint32_t nr_regions;
	struct region_shape_t region_shapes[MAX_REGIONS_PER_DPU];
};

struct dpu_work_t regions_processing[MAX_RANKS][MAX_DPUS_PER_RANK];

// No verification is done on the use of those buffers.
// If TOTAL_READS and TOTAL_HAPS are too small, data may get corrupted without any warning.
uint64_t reads_len_buffer       [TOTAL_READS];
char     reads_array_buffer     [TOTAL_READS * MAX_READ_LENGTH];
uint32_t priors                 [TOTAL_READS * MAX_READ_LENGTH];
int32_t  match_to_indel_buffer  [TOTAL_READS * MAX_READ_LENGTH];
uint64_t haplotypes_len_buffer  [TOTAL_HAPS];
uint64_t haplotypes_val_buffer  [TOTAL_HAPS];
char     haplotypes_array_buffer[TOTAL_HAPS * MAX_HAPLOTYPE_LENGTH];

struct dpu_regions_t dpu_regions_buffer[DPU_INPUT_BUFFER_SIZE];
struct queue_t       dpu_regions_queue;

struct dpu_results_t dpu_results_buffer[DPU_OUTPUT_BUFFER_SIZE];
struct queue_t       dpu_results_queue;
