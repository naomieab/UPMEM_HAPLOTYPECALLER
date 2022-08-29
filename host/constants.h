//#define DECIMAL_BITS 12
//#define ONE  512
#define BITS 12
#define DECIMAL_BITS 9
#define DECIMAL_MASK ( -1 >>> (BITS - DECIMAL_BITS) )
#define FULL_MASK (-1)
#define MAX_VALUE INT_MAX
#define MIN_VALUE INT_MIN
#define LUT_SIZE 1202 // should be 4096 (equals to 2^BITS) but in reality values bigger than 1202 are pruned to zero
#define ONE 512//(int)1 << DECIMAL_BITS

/* BITS_MASK is a mask composed of ones in place of the bits and zeros elsewhere
	for example for a word with 2 bits we'll have 30 zeros and them 2 ones
   UNBITS_MASK is the not of BITS_MASK
   The two masks must correspond to the bits number
*/
#define BITS_MASK 4095 // 0b0000 0000 0000 0000 00000 1111 1111 1111
#define UNBITS_MASK -4096 // 0b1111 1111 1111 1111 1111 0000 0000 0000
#define LIMIT 2048 //tranfer size limit for MRAM_READ


#define MAX_READ_LENGTH 80 //must be a multiple of 8 (because for the biggest read we bring READ_SIZE*sizeof(char))
#define MAX_READ_NUM 2500//1432
#define MAX_HAPLOTYPE_LENGTH 380//200 //MUST BE MULTIPLE OF 4
#define MAX_HAPLOTYPE_NUM 24//64


#define TOTAL_REGIONS 132108//167518//Total number of regions in the run
#define NR_REGIONS 2546 //Total number of dpus available (we send by chunk to dpus)

//sum of all the reads/haplotypes in chunks of 2546 regions ( + maximum number of reads/haplotypes per region)
#define TOTAL_READS 1165000//44853 
#define TOTAL_HAPS 62000//8000 

//dedicated to an array containing the offsets of the data
#define OFFSET_SIZE 5
#define READS_LEN_ARRAY 0
#define HAPLOTYPES_LEN_VAL_ARRAY 1
#define READS_ARR 2
#define HAPS_ARR 3
#define PRIOR_ARR 4

//void createLut(int[] lut);