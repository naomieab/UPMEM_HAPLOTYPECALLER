CC=gcc
DPU_DIR := dpu
HOST_DIR := host
BUILDDIR ?= build
NR_TASKLETS ?= 12
NR_DPUS ?= 1

HOST_TARGET := ${BUILDDIR}/haplotype_host
DPU_TARGET := ${BUILDDIR}/haplotype_dpu

HOSTFIXED_DIR := fixedSize
DPUFIXED_DIR := fixedSize/dpu/

HOSTF_SOURCES := $(wildcard ${HOSTFIXED_DIR}/*.c)
DPUF_SOURCES := $(wildcard ${DPUFIXED_DIR}/*.c)
HOSTF_TARGET := ${BUILDDIR}/haplotype_hostFixedSize
DPUF_TARGET := ${BUILDDIR}/haplotype_dpuFixedSize


HOST_SOURCES := $(wildcard ${HOST_DIR}/*.h ${HOST_DIR}/*.c)
DPU_SOURCES := $(wildcard ${DPU_DIR}/*.c)
DPU_DEPENDENCIES := $(wildcard ${DPU_DIR}/*.h ${DPU_DIR}/*.c)

.PHONY: all clean test_c


HOST_FLAGS := -std=c11 -O3 -lm `dpu-pkg-config --cflags --libs dpu` -DNR_TASKLETS=${NR_TASKLETS} -DNR_DPUS=${NR_DPUS}
DPU_FLAGS := -DNR_TASKLETS=${NR_TASKLETS}

all: ${HOST_TARGET} ${DPU_TARGET}


fixedSize: ${HOSTF_TARGET} ${DPUF_TARGET}
	./${HOSTF_TARGET}

${HOSTF_TARGET}: ${HOSTF_SOURCES}
	$(CC) -o $@ ${HOSTF_SOURCES} ${HOST_FLAGS}

${DPUF_TARGET}: ${DPUF_SOURCES} 
	dpu-upmem-dpurte-clang ${DPU_FLAGS} -o $@ ${DPUF_SOURCES}




${HOST_TARGET}: ${HOST_SOURCES} 
	$(CC) -o $@ ${HOST_SOURCES} ${HOST_FLAGS}

${DPU_TARGET}: ${DPU_DEPENDENCIES} 
	dpu-upmem-dpurte-clang ${DPU_FLAGS} -o $@ ${DPU_SOURCES}

clean:
	$(RM) -r $(BUILDDIR)

test_valgrind: ${HOST_TARGET} ${DPU_TARGET}
	valgrind -s --leak-check=full ./${HOST_TARGET} ${INPUT} ${OUTPUT} ${PERF}

test: ${HOST_TARGET} ${DPU_TARGET}
	./${HOST_TARGET} ${INPUT} ${OUTPUT} ${PERF}
