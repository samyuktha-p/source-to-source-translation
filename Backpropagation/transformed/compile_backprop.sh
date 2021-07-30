#!/bin/bash 
# /home/samyuktha/Desktop/compiler/cetus-1.4.4/build.sh bin
cetus -opencl-dpu -outdir=. -profitable-omp=0 ./original/backprop_ocl.c ./original/backprop_kernel_1.cl ./original/backprop_kernel_2.cl
#dpu-upmem-dpurte-clang -DNR_TASKLETS=${1} -O2 -o bpnn_adjust_weights_ocl bpnn_adjust_weights_ocl.c
#dpu-upmem-dpurte-clang -DNR_TASKLETS=${1} -O2 -o bpnn_layerforward_ocl bpnn_layerforward_ocl.c
#rm backprop_host
#make
