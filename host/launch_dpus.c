#include "buffers.h"
#include "circular_queue.h"
#include "constants.h"
#include "log.h"
#include <dpu.h>
#include <stdint.h>

#ifndef DPU_BINARY
#define DPU_BINARY "./build/haplotype_dpu"
#endif

static struct region_shape_t empty_region_shape = {0,0,0};
static int64_t dummy_likelihoods[MAX_READ_NUM][MAX_HAPLOTYPE_NUM];

static pthread_mutex_t populating_mutex = PTHREAD_MUTEX_INITIALIZER;

dpu_error_t run_rank(struct dpu_set_t set, uint32_t rank_id, void* arg) {
	uint32_t each_dpu;
	struct dpu_set_t dpu;
	int dpu_buffer_indices[MAX_DPUS_PER_RANK];
	// Reserve all needed regions for this rank.
	bool all_dpus_inactive = true;
	bool queue_closed = false;
	pthread_mutex_lock(&populating_mutex);
	DPU_FOREACH(set, dpu, each_dpu) {
		int region_id;
		if (queue_closed) {
			region_id = -1;
		} else {
			LOG_DEBUG("taking from queue\n");
			region_id = queue_take(&dpu_regions_queue);
		}
		dpu_buffer_indices[each_dpu] = region_id;
		if (region_id < 0) {
			queue_closed = true;
		} else {
			regions_processing[rank_id][each_dpu].first_region_index = dpu_regions_buffer[region_id].first_region_index;
			regions_processing[rank_id][each_dpu].nr_regions = dpu_regions_buffer[region_id].nr_regions;
			memcpy(regions_processing[rank_id][each_dpu].region_shapes, dpu_regions_buffer[region_id].region_shapes, sizeof(struct region_shape_t) * MAX_REGIONS_PER_DPU);
			all_dpus_inactive = false;
		}
	}
	if (all_dpus_inactive) {
		queue_close_writer(&dpu_results_queue, 100);// FIXME: get actual max number of readers for that queue.
		return DPU_OK;
	}

	DPU_FOREACH(set, dpu, each_dpu) {
		if (dpu_buffer_indices[each_dpu] < 0) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, &priors));
		} else {
			DPU_ASSERT(dpu_prepare_xfer(dpu, &dpu_regions_buffer[dpu_buffer_indices[each_dpu]].nr_haplotypes));
		}
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "nr_haplotypes", 0, sizeof(uint64_t), DPU_XFER_DEFAULT));

	DPU_FOREACH(set, dpu, each_dpu) {
		if (dpu_buffer_indices[each_dpu] < 0) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, &priors));
		} else {
			DPU_ASSERT(dpu_prepare_xfer(dpu, &dpu_regions_buffer[dpu_buffer_indices[each_dpu]].nr_reads));
		}
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "nr_reads", 0, sizeof(uint64_t), DPU_XFER_DEFAULT));

	DPU_FOREACH(set, dpu, each_dpu) {
		if (dpu_buffer_indices[each_dpu] < 0) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, &priors));
		} else {
			DPU_ASSERT(dpu_prepare_xfer(dpu, &dpu_regions_buffer[dpu_buffer_indices[each_dpu]].read_region_starts));
		}
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "read_region_starts", 0, (MAX_REGIONS_PER_DPU+1) * sizeof(uint32_t), DPU_XFER_DEFAULT));

	DPU_FOREACH(set, dpu, each_dpu) {
		if (dpu_buffer_indices[each_dpu] < 0) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, &priors));
		} else {
			DPU_ASSERT(dpu_prepare_xfer(dpu, &dpu_regions_buffer[dpu_buffer_indices[each_dpu]].haplotype_region_starts));
		}
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "haplotype_region_starts", 0, (MAX_REGIONS_PER_DPU+1) * sizeof(uint32_t), DPU_XFER_DEFAULT));

	DPU_FOREACH(set, dpu, each_dpu) {
		if (dpu_buffer_indices[each_dpu] < 0) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, &priors));
		} else {
			DPU_ASSERT(dpu_prepare_xfer(dpu, dpu_regions_buffer[dpu_buffer_indices[each_dpu]].reads_len));
		}
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "mram_reads_len", 0, MAX_READ_NUM * sizeof(uint64_t), DPU_XFER_DEFAULT));

	DPU_FOREACH(set, dpu, each_dpu) {
		if (dpu_buffer_indices[each_dpu] < 0) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, &priors));
		} else {
			DPU_ASSERT(dpu_prepare_xfer(dpu, dpu_regions_buffer[dpu_buffer_indices[each_dpu]].haplotypes_len));
		}
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "mram_haplotypes_len", 0, MAX_HAPLOTYPE_NUM * sizeof(uint64_t), DPU_XFER_DEFAULT));

	DPU_FOREACH(set, dpu, each_dpu) {
		if (dpu_buffer_indices[each_dpu] < 0) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, &priors));
		} else {
			DPU_ASSERT(dpu_prepare_xfer(dpu, dpu_regions_buffer[dpu_buffer_indices[each_dpu]].reads_array));
		}
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "mram_reads_array", 0, MAX_READ_NUM * MAX_READ_LENGTH * sizeof(char), DPU_XFER_DEFAULT));

	DPU_FOREACH(set, dpu, each_dpu) {
		if (dpu_buffer_indices[each_dpu] < 0) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, &priors));
		} else {
			DPU_ASSERT(dpu_prepare_xfer(dpu, dpu_regions_buffer[dpu_buffer_indices[each_dpu]].match_to_indel));
		}
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "mram_matchToIndelArray", 0, MAX_READ_NUM * MAX_READ_LENGTH * sizeof(int32_t), DPU_XFER_DEFAULT));

	DPU_FOREACH(set, dpu, each_dpu) {
		if (dpu_buffer_indices[each_dpu] < 0) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, &priors));
		} else {
			DPU_ASSERT(dpu_prepare_xfer(dpu, dpu_regions_buffer[dpu_buffer_indices[each_dpu]].priors));
		}
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "mram_priors", 0, 2 * MAX_READ_NUM * MAX_READ_LENGTH * sizeof(int32_t), DPU_XFER_DEFAULT));

	DPU_FOREACH(set, dpu, each_dpu) {
		if (dpu_buffer_indices[each_dpu] < 0) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, &priors));
		} else {
			DPU_ASSERT(dpu_prepare_xfer(dpu, dpu_regions_buffer[dpu_buffer_indices[each_dpu]].haplotypes_array));
		}
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "mram_haplotypes_array", 0, MAX_HAPLOTYPE_NUM * MAX_HAPLOTYPE_LENGTH * sizeof(char), DPU_XFER_DEFAULT));

	DPU_FOREACH(set, dpu, each_dpu) {
		if (dpu_buffer_indices[each_dpu] < 0) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, &priors));
		} else {
			DPU_ASSERT(dpu_prepare_xfer(dpu, dpu_regions_buffer[dpu_buffer_indices[each_dpu]].haplotypes_val));
		}
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_TO_DPU, "mram_haplotypes_val", 0, MAX_HAPLOTYPE_NUM * sizeof(uint64_t), DPU_XFER_DEFAULT));


	// Ensure all data transfers are done before releasing it.
	//LOG_DEBUG("dpu sync\n");
	//DPU_ASSERT(dpu_sync(set));
	// Release all regions sent.
	DPU_FOREACH(set, dpu, each_dpu) {
		queue_release(&dpu_regions_queue, dpu_buffer_indices[each_dpu]);
	}
	pthread_mutex_unlock(&populating_mutex);
	
	LOG_INFO("launching rank %d\n", rank_id);
	DPU_ASSERT(dpu_launch(set, DPU_SYNCHRONOUS));

	//Get results back
	LOG_INFO("fetching results for rank %d\n", rank_id);

	DPU_FOREACH(set, dpu, each_dpu) {
		if (dpu_buffer_indices[each_dpu] >= 0) {
			int result_id = queue_put(&dpu_results_queue);
			dpu_buffer_indices[each_dpu] = result_id;
			DPU_ASSERT(dpu_prepare_xfer(dpu, dpu_results_buffer[dpu_buffer_indices[each_dpu]].likelihoods));
		} else {
			DPU_ASSERT(dpu_prepare_xfer(dpu, dummy_likelihoods));
		}
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_FROM_DPU, "likelihoods", 0, MAX_HAPLOTYPE_NUM * MAX_READ_NUM * sizeof(int64_t), DPU_XFER_DEFAULT));

	uint64_t dummy_nr_cycles;
	DPU_FOREACH(set, dpu, each_dpu) {
		if (dpu_buffer_indices[each_dpu] < 0) {
			DPU_ASSERT(dpu_prepare_xfer(dpu, &dummy_nr_cycles));
		} else {
			DPU_ASSERT(dpu_prepare_xfer(dpu, &(dpu_results_buffer[dpu_buffer_indices[each_dpu]].nr_cycles)));
		}
	}
	DPU_ASSERT(dpu_push_xfer(set, DPU_XFER_FROM_DPU, "nb_cycles", 0, sizeof(uint64_t), DPU_XFER_DEFAULT));

	DPU_FOREACH(set, dpu, each_dpu) {
		if (dpu_buffer_indices[each_dpu] >= 0) {
			dpu_results_buffer[dpu_buffer_indices[each_dpu]].first_region_index = regions_processing[rank_id][each_dpu].first_region_index;
			dpu_results_buffer[dpu_buffer_indices[each_dpu]].nr_regions = regions_processing[rank_id][each_dpu].nr_regions;
			memcpy(dpu_results_buffer[dpu_buffer_indices[each_dpu]].region_shapes, regions_processing[rank_id][each_dpu].region_shapes, sizeof(struct region_shape_t) * MAX_REGIONS_PER_DPU);
		}
	}

	DPU_ASSERT(dpu_sync(set));
	DPU_FOREACH(set, dpu, each_dpu) {
		if (dpu_buffer_indices[each_dpu] >= 0) {
			queue_make_available(&dpu_results_queue, dpu_buffer_indices[each_dpu]);
		}
	}

	
	if (queue_closed) {
		queue_close_writer(&dpu_results_queue, 100);// FIXME: get actual max number of readers for that queue.
		return DPU_OK;
	}
	// Callback itself
	DPU_ASSERT(dpu_callback(set, run_rank, NULL, DPU_CALLBACK_ASYNC));
	return DPU_OK;
}

static struct dpu_set_t set_all_dpus;

void launch_all_ranks() {
	int nr_ranks;
	int nr_dpus;

	DPU_ASSERT(dpu_alloc(DPU_ALLOCATE_ALL, NULL, &set_all_dpus));
	DPU_ASSERT(dpu_load(set_all_dpus, DPU_BINARY, NULL));
	DPU_ASSERT(dpu_get_nr_ranks(set_all_dpus, &nr_ranks));
	DPU_ASSERT(dpu_get_nr_dpus(set_all_dpus, &nr_dpus));
	LOG_INFO("Nr ranks =%d, and Nr dpus = %d\n", nr_ranks, nr_dpus);

	for (int i=0; i<nr_ranks; i++) {
		queue_open_writer(&dpu_results_queue);
	}

	DPU_ASSERT(dpu_callback(set_all_dpus, run_rank, NULL, DPU_CALLBACK_ASYNC));
}

void free_dpus() {
	DPU_ASSERT(dpu_free(set_all_dpus));
}
