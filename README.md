# UPMEM_HAPLOTYPECALLER
An implementation of the GATK HaplotypeCaller algorithm using UPMEM chip (in memory processing)

The first implementation used the original GATK HC algorithm, but running this algorithm on the UPMEM DPU 
makes it slower since it involves a lot of floating point calculations as well as multiplications (which 
are costly for DPUs). 

With this observation we decided to change the algorithm to use fixed precision integer. Such a change 
doesn't impact much the accuracy/precision of the results on a whole genome sequencing. 

The second decision is to perform the calculations in the log domain. 
This decision has 2 consequences: first and most important, multiplications become additions, and DPUs have
better performances with additions; second, this reduces the range of values making possible the use of a 
LUT (of moderate size) who's entry are the bits of the fixed point integer. 

The implementation is divided between the DPU code (under dpu) and the host code (under host)

A makefile is also added (NR_TASKLETS must be a multiple of 2)

Two scripts are added to the repository: 
-regionStats.c computes certain statistics on an input regions files such as the number of regions, the maximum number of reads/haplotypes, the maximum length of reads/haplotypes...
It also print to the output file provided (csv file) for each region the computation complexity of the region, as well as the number of reads and haplotypes in the region. 
Command line to run:
./regionStats input_file.csv output_file.csv 
The input file is a file containing the regions details (haplotype number / haplotype list / read number / read+quality list)
Once this script is ran you must open the output file (excel) and sort according to column B (this new file will serve as data to the next script)

-sortRegions.c creates a new regions file that is sorted and splits last regions. 
In order to run it you must first complete the sortRegions.h: 
  -update NR_REGIONS to the relevant number of regions 
  -update TO_SEPARATE to the number of regions you wish to split
  -update the array named sorted as following: the array must contain the A column of the output file of the previous script (ie, the order of the  regions)
  (copy this column to notepad++ for example and insert a comma at each line (in notepad++ you can search $ and replace by ,)
Command line to run:
./sortRegions input_file.csv output_file.csv
input_file is the non sorted file of regions and output_file is the sorted/split file of regions


Once these steps have been completed you can run the HaplotypeCaller algorithm on DPUs with the following line:
make test input_file.csv output_file.csv performance_file.csv
