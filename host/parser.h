#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <dpu_types.h>

#define BUFFER_SIZE 4098

void read_and_allocate(FILE* file);

void add_haplotype(FILE* file, int region, int index);

void add_read(FILE* file, int region, int index);

FILE* read_data(char* filename);

void free_mem(FILE* file);