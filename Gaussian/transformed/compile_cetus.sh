cetus -opencl-dpu -outdir=. -profitable-omp=0 ./original/gaussianElim.c ./original/clutils.c ./original/utils.c ./original/gaussianElim_kernels_1.cl ./original/gaussianElim_kernels_2.cl
#dpu-upmem-dpurte-clang -DNR_TASKLETS=${1} -O2 -o Fan1 Fan1.c
#dpu-upmem-dpurte-clang -DNR_TASKLETS=${1} -O2 -o Fan2 Fan2.c
#rm gaussian_host
#make KERNEL_DIM="-DRD_WG_SIZE_0=${2} -DRD_WG_SIZE_1_0=${2} -DRD_WG_SIZE_1_1=${2}"
