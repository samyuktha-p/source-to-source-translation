#!/bin/bash
for g in 1 2 3 4 5 6
do
	t=$g
	while [ $t -le 22 ]
	do
		echo "$g - $t"
		dpu-upmem-dpurte-clang -DNR_TASKLETS=${t} -O2 -o Fan1 Fan1.c
		dpu-upmem-dpurte-clang -DNR_TASKLETS=${t} -O2 -o Fan2 Fan2.c
		rm gaussian_host
		make KERNEL_DIM="-DRD_WG_SIZE_0=${1} -DRD_WG_SIZE_1_0=${1} -DRD_WG_SIZE_1_1=${1}"
		./gaussian_host -f ../../data/gaussian/matrix208.txt -n $(($t)) -g $(($g))
		t=$(($t+$g))
	done
done
