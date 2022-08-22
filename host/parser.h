#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <dpu_types.h>
#include <math.h>
#include <stdbool.h>

#define BUFFER_SIZE 4098


void add_haplotype(FILE* file, int hap_idx, int index);

void add_read(FILE* file, int read_idx, int index);

FILE* read_data(FILE* file, int nr_dpus);

void free_mem(FILE* file);
