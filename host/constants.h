
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
#define UNBITS_MASK 4294963200 // 0b1111 1111 1111 1111 0000 0000 0000

#define MAX_READ_LENGTH 120
#define MAX_READ_NUM 100
#define MAX_HAPLOTYPE_LENGTH 80
#define MAX_HAPLOTYPE_NUM 16

#define NR_REGIONS 64

//void createLut(int[] lut);