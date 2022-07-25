#include <dpu_types.h>
#include <stdlib.h>
#include <math.h>
#include <dpu.h>
#include <assert.h>

void free_prior(int** array, int region);
void create_prior();
void populate_mram(struct dpu_set_t set, uint32_t nr_dpus, int region);
