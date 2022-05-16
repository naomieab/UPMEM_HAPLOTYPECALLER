#include <mram.h>
#include <defs.h>
#include <stdio.h>
#include <perfcounter.h>
#include <alloc.h>
#include <barrier.h>
#include <sysdef.h>
#include <assert.h>
#include <string.h>

#define DEFAULT_INSERTION_DELETION_QUAL 45
#define gcpHMM 10
#define QUALITY_OFFSET 33

#define matchToMatch 0
#define indelToMatch 1
#define matchToInsertion 2
#define insertionToInsertion 3
#define matchToDeletion 4
#define deletionToDeletion 5
#define TRANS_PROB_ARRAY_LENGTH 6

#define NR_HAPLOTYPE 1
#define HAPLOTYPE_LEN 10
#define MATRIX_LINES 2 //number of lines we need to store for each matrix
#define READ_LEN 120
#define NR_READ 16

#define MAX_READS 250
#define MAX_HAPLOTYPES 100




