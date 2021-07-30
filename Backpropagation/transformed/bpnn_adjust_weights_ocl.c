#include<stdio.h>
#include <alloc.h>
#include <mram.h>
#include <stddef.h>
#include <defs.h>
#include <mutex.h>
#include <barrier.h>
#include <perfcounter.h>
#include <stdint.h>

__host __dma_aligned int64_t N_TASKLETS = NR_TASKLETS;
#define DPU_MRAM_SIZE 10108864
MUTEX_INIT(my_mutex);
BARRIER_INIT(my_barrier, NR_TASKLETS);

__host __dma_aligned int64_t cycles = 0;
__host __dma_aligned int64_t INPUT = 0;
__host __dma_aligned int64_t GRID_LIMIT = 1;

//__host __dma_aligned int64_t barrier_count[NR_TASKLETS];
//__host __dma_aligned int64_t T_TD_START[NR_TASKLETS];
//__host __dma_aligned int64_t T_TD_END[NR_TASKLETS];
//__host __dma_aligned int64_t N_WG_ID[NR_TASKLETS];
__host __dma_aligned int64_t WORK_DONE[NR_TASKLETS];

__host __dma_aligned int64_t MULTI_WGS = 0;

__host __dma_aligned int64_t _bd[] = {0, 0, 0};
__dma_aligned int64_t start_gd[] = {0, 0, 0};
__host __dma_aligned int64_t _gs[] = {1, 1, 1};
__host __dma_aligned int64_t _bs[] = {1, 1, 1};

__host __dma_aligned int64_t dpu_gs[] = {1, 0, 0};

int B_THREADS = 0;
int T_THREADS = 0;

__host __dma_aligned int64_t PARTITION_DIM_GRID = 1;
__host __dma_aligned int64_t PARTITION_DIM_WG = 0;

int NTASKLET_WG[2] = {0 , 0};
int NTASK_TASKLET[2][2] = {{0, 0}, {0, 0}};

typedef struct {
	int start;
	int end;
} Point;
__host __dma_aligned int64_t t_dimension = 2;

BARRIER_INIT(barrier_0, 16);
// BARRIER_INIT(barrier_1, 1);
// BARRIER_INIT(barrier_2, 1);
// BARRIER_INIT(barrier_3, 1);
// BARRIER_INIT(barrier_4, 1);
// BARRIER_INIT(barrier_5, 1);
// BARRIER_INIT(barrier_6, 1);
// BARRIER_INIT(barrier_7, 1);
// BARRIER_INIT(barrier_8, 1);
// BARRIER_INIT(barrier_9, 1);
// BARRIER_INIT(barrier_10, 1);
// BARRIER_INIT(barrier_11, 1);
// BARRIER_INIT(barrier_12, 1);
// BARRIER_INIT(barrier_13, 1);
// BARRIER_INIT(barrier_14, 1);
// BARRIER_INIT(barrier_15, 1);

// barrier_t* barrier[] = {&barrier_0, &barrier_1, &barrier_2, &barrier_3, &barrier_4, &barrier_5, &barrier_6, &barrier_7, &barrier_8, &barrier_9, &barrier_10, &barrier_11, &barrier_12, &barrier_13, &barrier_14, &barrier_15};
barrier_t* barrier[] = {&barrier_0};
#define N_MULTI_WGS 1

__host __dma_aligned long long oldw_map[2];
__host __dma_aligned long long oldw_stride[2];
__host __dma_aligned long long oldw_length[2];
__host __dma_aligned long long oldw_flag[2];
__host __dma_aligned long long oldw_cpu_offset[2];
__host __dma_aligned long long oldw_offset[2];
__host __dma_aligned long long w_map[2];
__host __dma_aligned long long w_stride[2];
__host __dma_aligned long long w_length[2];
__host __dma_aligned long long w_flag[2];
__host __dma_aligned long long w_cpu_offset[2];
__host __dma_aligned long long w_offset[2];
__host __dma_aligned long long ly_map[1];
__host __dma_aligned long long ly_stride[1];
__host __dma_aligned long long ly_length[1];
__host __dma_aligned long long ly_flag[1];
__host __dma_aligned long long ly_cpu_offset[1];
__host __dma_aligned long long ly_offset[1];
__host __dma_aligned long long delta_map[1];
__host __dma_aligned long long delta_stride[1];
__host __dma_aligned long long delta_length[1];
__host __dma_aligned long long delta_flag[1];
__host __dma_aligned long long delta_cpu_offset[1];
__host __dma_aligned long long delta_offset[1];
int max(int a, int b)
{
	return (a<b) ?b :a;
}

int min(int a, int b)
{
	return (a<b) ?a :b;
}

void initialize_BT()
{
	for(int dim=0; dim<3; dim++) {
	 start_gd[dim] = _bd[dim] * _bs[dim];
}

}

void divide_work_group(int category)
{
if(PARTITION_DIM_WG==-1)
	return;

int n_TASKS = _bs[PARTITION_DIM_WG];
int n_TASKLETS_WG = NTASKLET_WG[0] + (1-category);

NTASK_TASKLET[category][0] = n_TASKS/n_TASKLETS_WG;
NTASK_TASKLET[category][1] = n_TASKS%n_TASKLETS_WG;
}

void divide_grid()
{
if(PARTITION_DIM_GRID==-1) {
	NTASKLET_WG[0] = NR_TASKLETS;
}
else {
	NTASKLET_WG[0] = NR_TASKLETS/dpu_gs[PARTITION_DIM_GRID];
	NTASKLET_WG[1] = NR_TASKLETS%dpu_gs[PARTITION_DIM_GRID];
}

if(NTASKLET_WG[1]!=0) {
	divide_work_group(0);
}
	divide_work_group(1);
}

int get_WORK_GROUP_ID(int TASKLET_ID, int * t_bd, int * TASKLET_ID_WG)
{
if(PARTITION_DIM_GRID==-1)
	return 1;

int TBREAK_GRID = NR_TASKLETS - NTASKLET_WG[1];

if(TASKLET_ID < TBREAK_GRID) {
	t_bd[PARTITION_DIM_GRID] = TASKLET_ID/NTASKLET_WG[0];
	*TASKLET_ID_WG = TASKLET_ID%NTASKLET_WG[0];

if( t_bd[0] < NTASKLET_WG[1] )
	return 0;
else
	return 1;
}

else {
	t_bd[PARTITION_DIM_GRID] = TASKLET_ID - TBREAK_GRID;
	*TASKLET_ID_WG = NTASKLET_WG[0];
	return 0;
}
}

void get_start_end(int TASKLET_ID_WG, int category, Point * t_td_se)
{
if(PARTITION_DIM_WG == -1)
	return;

int n_TASKLETS_WG = NTASKLET_WG[0] + (1-category);

t_td_se[PARTITION_DIM_WG].start = NTASK_TASKLET[category][0]*TASKLET_ID_WG;

if(TASKLET_ID_WG < NTASK_TASKLET[category][1]) {
	t_td_se[PARTITION_DIM_WG].start += TASKLET_ID_WG;
	t_td_se[PARTITION_DIM_WG].end = 1;
}
else {
	t_td_se[PARTITION_DIM_WG].start += NTASK_TASKLET[category][1];
	t_td_se[PARTITION_DIM_WG].end = 0;
}
t_td_se[PARTITION_DIM_WG].end += t_td_se[PARTITION_DIM_WG].start + NTASK_TASKLET[category][0];
}

int get_relative_global_id(int dim, int WG_ID, int _td_dim)
{
int rem_threads = 0;
	if(dim==PARTITION_DIM_GRID)
		 rem_threads += WG_ID*_bs[dim];
rem_threads += _td_dim;

	 return rem_threads;

}

int get_global_id(int dim, int WG_ID, int _td_dim)
{
return start_gd[dim]+get_relative_global_id(dim, WG_ID, _td_dim);

}

void print_dpu_mat( char * arr_name,  __mram_ptr float * arr,  int size)
{
printf("%s:: ", arr_name);
for(int i=0; i<size; i++) {
printf("%f ", arr[i]);
}
printf("\n");
}

/*
Copyright (C) 1991-2018 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it andor
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <http:www.gnu.org/licenses/>. 
*/
/*
This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it. 
*/
/*
glibc's intent is to support the IEC 559 math functionality, real
   and complex.  If the GCC (4.9 and later) predefined macros
   specifying compiler intent are available, use them to determine
   whether the overall intent is to support these features; otherwise,
   presume an older compiler has intent to support these features and
   define these macros by default. 
*/
/*
wchar_t uses Unicode 10.0.0.  Version 10.0 of the Unicode Standard is
   synchronized with ISOIEC 10646:2017, fifth edition, plus
   the following additions from Amendment 1 to the fifth edition:
   - 56 emoji characters
   - 285 hentaigana
   - 3 additional Zanabazar Square characters
*/
/* We do not support C11 <threads.h>.  */
__host __dma_aligned __mram_ptr float * delta;
__host __dma_aligned int hid;
__host __dma_aligned __mram_ptr float * ly;
__host __dma_aligned int in;
__host __dma_aligned __mram_ptr float * w;
__host __dma_aligned __mram_ptr float * oldw;
#define SIZE_DEFAULT 100

__host __dma_aligned int64_t delta_start = SIZE_DEFAULT;
__host __dma_aligned int64_t ly_start = SIZE_DEFAULT;
__host __dma_aligned int64_t w_start = SIZE_DEFAULT;
__host __dma_aligned int64_t oldw_start = SIZE_DEFAULT;

void bpnn_adjust_weights_ocl(int INPUT, int TASKLET_ID_WG, int * t_bd, Point * t_td_se)
{
int tid = me();
int ss_tx = 0, ss_ty = 0, ss_tz = 0;
int by = t_bd[1];
int tx = tx;
int ty = ty;
/* int index =  ( hid + 1 ) HEIGHT * by + ( hid + 1 ) * ty + tx + 1 + ( hid + 1 ) ;   */
/* int index_y = HEIGHT by + ty + 1; */
/* int index_x = tx + 1; */
if(tid == 0) {
delta = DPU_MRAM_HEAP_POINTER + delta_start; 
ly = DPU_MRAM_HEAP_POINTER + ly_start; 
w = DPU_MRAM_HEAP_POINTER + w_start; 
oldw = DPU_MRAM_HEAP_POINTER + oldw_start; 
}
if(TASKLET_ID_WG == 0) {
}
barrier_wait(&my_barrier);
if(_bd[PARTITION_DIM_GRID]+t_bd[PARTITION_DIM_GRID]>=GRID_LIMIT) {
// printf("###SKIP### %d %lld %lld\n", me(), _bd[PARTITION_DIM_GRID]+t_bd[PARTITION_DIM_GRID], GRID_LIMIT);
	int BARRIER_COUNT = 1;
	for(int ss_b=0; ss_b<BARRIER_COUNT; ss_b++)
		barrier_wait(barrier[t_bd[PARTITION_DIM_GRID]]);
	return;}
;
for (ty=(t_td_se[1].start); ty<(t_td_se[1].end); ty ++ )
{
for (tx=(t_td_se[0].start); tx<(t_td_se[0].end); tx ++ )
{
w[((((((hid+1)*16)*by)+((hid+1)*ty))+tx)+0)+0]+=(((0.3F*delta[tx+0])*ly[((16*by)+ty)+0])+(0.3F*oldw[((((((hid+1)*16)*by)+((hid+1)*ty))+tx)+0)+0]));
}
}
for (ty=(t_td_se[1].start); ty<(t_td_se[1].end); ty ++ )
{
for (tx=(t_td_se[0].start); tx<(t_td_se[0].end); tx ++ )
{
oldw[((((((hid+1)*16)*by)+((hid+1)*ty))+tx)+0)+0]=(((0.3F*delta[tx+0])*ly[((16*by)+ty)+0])+(0.3F*oldw[((((((hid+1)*16)*by)+((hid+1)*ty))+tx)+0)+0]));
}
}
barrier_wait(barrier[t_bd[PARTITION_DIM_GRID]]);
//++barrier_count[tid];

for (ty=(t_td_se[1].start); ty<(t_td_se[1].end); ty ++ )
{
if ((ty==0)&&(by==0))
{
for (tx=(t_td_se[0].start); tx<(t_td_se[0].end); tx ++ )
{
w[(tx+0)+w_offset[1]]+=((0.3F*delta[tx+0])+(0.3F*oldw[(tx+0)+oldw_offset[1]]));
}
for (tx=(t_td_se[0].start); tx<(t_td_se[0].end); tx ++ )
{
oldw[(tx+0)+oldw_offset[1]]=((0.3F*delta[tx+0])+(0.3F*oldw[(tx+0)+oldw_offset[1]]));
}
}
}
return ;
}

int main()
{
if(INPUT==-1) {
cycles = (uint32_t)1;
return 0;
}

if(INPUT==0) {
if(me()==0) {
	divide_grid();
	perfcounter_config(COUNT_CYCLES, true);
++INPUT;
	//printf("%d - NTASKLET_WG = %d %d\n", me(), NTASKLET_WG[0], NTASKLET_WG[1]);
	//printf("%d - NTASK_TASKLET = %d-%d, %d-%d\n", me(), NTASK_TASKLET[0][0], NTASK_TASKLET[0][1], NTASK_TASKLET[1][0], NTASK_TASKLET[1][1]);
}
barrier_wait(&my_barrier);

return 0;
}
else if(INPUT==1) {
	if(me()==0)
		initialize_BT();
	barrier_wait(&my_barrier);

}
int TASKLET_ID_WG = me();

int t_bd[3] = {_bd[0], _bd[1], _bd[2]};
Point t_td_se[3] = {{0, _bs[0]}, {0, _bs[1]}, {0, _bs[2]}};

int category = get_WORK_GROUP_ID(me(), t_bd, &TASKLET_ID_WG);

//printf("%d - TASKLET_ID_WG = %d\n", me(), TASKLET_ID_WG);
//printf("%d - t_bd = %d %d %d\n", me(), t_bd[0], t_bd[1], t_bd[2]);
get_start_end(TASKLET_ID_WG, category, t_td_se);
//printf("%d - t_td_se = %d-%d, %d-%d, %d-%d\n", me(), t_td_se[0].start, t_td_se[0].end, t_td_se[1].start, t_td_se[1].end, t_td_se[2].start, t_td_se[2].end);
 //T_TD_START[me()] = t_td_se[PARTITION_DIM_WG].start;
//T_TD_END[me()] = t_td_se[PARTITION_DIM_WG].end;
//N_WG_ID[me()] = TASKLET_ID_WG;
bpnn_adjust_weights_ocl(0, TASKLET_ID_WG, t_bd, t_td_se);

if(me()==0) {
	cycles = (uint32_t)perfcounter_get();
	printf("cycles = %lld\n", cycles);
}
barrier_wait(&my_barrier);

WORK_DONE[me()] = 1;

return 0;
}

