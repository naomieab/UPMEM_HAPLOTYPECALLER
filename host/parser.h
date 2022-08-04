#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <dpu_types.h>
#include <math.h>
#include <stdbool.h>

#define BUFFER_SIZE 4098


void add_haplotype(FILE* file, int region, int index);

void add_read(FILE* file, int region, int index);

FILE* read_data(char* filename);

void free_mem(FILE* file);