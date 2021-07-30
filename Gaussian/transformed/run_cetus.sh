#!/bin/bash
dpu-upmem-dpurte-clang -DNR_TASKLETS=${1} -O2 -o Fan1 Fan1.c
dpu-upmem-dpurte-clang -DNR_TASKLETS=${1} -O2 -o Fan2 Fan2.c
rm gaussian_host
make KERNEL_DIM="-DRD_WG_SIZE_0=${2} -DRD_WG_SIZE_1_0=${2} -DRD_WG_SIZE_1_1=${2}"
