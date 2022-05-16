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

A makefile is also added
