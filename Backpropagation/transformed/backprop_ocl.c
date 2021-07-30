#include<dpu.h>
#include<inttypes.h>
#include<math.h>
#include<limits.h>
void note_down(int* p_val, int kernel_num);
int make8ByteAligned(int num, int e_size);
int ss_min(int a, int b);
int ss_max(int a, int b);
int* generateRectangle(long long *start, long long *end, long long *stride, long long *map, long long *flag, int fstride, int size);
int ss_unify(int arr[][4], long long * map, long long* flag, int size);
void combineRect(int arr[][4], long long *flag, int i, int j);
int ss_isIntersect(int arr[][4], int i, int j);
long long TOTAL_COUNT_KA=0, TOTAL_COUNT_KF=0;
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
/* includes, system */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/* #include <math.h> */
#include <sys/time.h>
#include "backprop.h"
/*
#ifdef NV NVIDIA

	#include <oclUtils.h>

#else 

	#include <CLcl.h>

#endif
*/
#include "open_cl.h"
/*  */
/* local variables */
static cl_context context;
static cl_command_queue cmd_queue;
static cl_device_type device_type;
static cl_device_id * device_list;
static cl_int num_devices;
static int initialize(int use_gpu)
{
	cl_int result;
	size_t size;
	/* create OpenCL context */
	cl_platform_id platform_id;
cl_context_properties ctxprop[] = {0, (cl_context_properties)platform_id, 0};
	int _ret_val_0;
	if (CL_SUCCESS!=0)
	{
		printf("ERROR: clGetPlatformIDs(1,*,0) failed\n");
		_ret_val_0=( - 1);
		return _ret_val_0;
	}
	device_type=(use_gpu ? 0 : 0);
	context=CL_SUCCESS;
	/* get the list of GPUs */
	size=sizeof(cl_device_id);
	result=CL_SUCCESS;
	num_devices=((int)(size/sizeof (cl_device_id)));
	printf("num_devices = %d\n", num_devices);
	/* device_list = new cl_device_id[num_devices]; */
	device_list=((cl_device_id * )malloc(sizeof (cl_device_id)*num_devices));
	result=CL_SUCCESS;
	/* create command queue for the first device */
	cmd_queue=CL_SUCCESS;
	_ret_val_0=0;
	return _ret_val_0;
}

static int shutdown()
{
	/* release resources */
	int _ret_val_0;
	if (cmd_queue)
	{
		CL_SUCCESS;
	}
	if (context)
	{
		CL_SUCCESS;
	}
	/* if( device_list ) delete[] device_list; */
	if (device_list)
	{
		free(device_list);
	}
	/* reset all variables */
	cmd_queue=0;
	context=0;
	device_list=0;
	num_devices=0;
	device_type=0;
	_ret_val_0=0;
	return _ret_val_0;
}

double gettime()
{
	struct timeval t;
	double _ret_val_0;
	gettimeofday( & t, (void * )0);
	_ret_val_0=(t.tv_sec+(t.tv_usec*1.0E-6));
	return _ret_val_0;
}

unsigned int num_threads = 0;
unsigned int num_blocks = 0;
/*  */
/* Program main */
/*  */
int main(int argc, char * * argv)
{
	int _ret_val_0;
	setup(argc, argv);
	return _ret_val_0;
}

int * generateRectangle(long long * start, long long * end, long long * stride, long long * map, long long * flag, int fstride, int size)
{
	int (*arr)[4] = (int (*)[4])malloc(size*sizeof(int)*4);
	for(int i=0; i<size; i++) {
		arr[i][0] = start[i]%fstride;
		arr[i][1] = start[i]/fstride;
		arr[i][2] = end[i]%fstride;
		arr[i][3] = end[i]/fstride;
		if(arr[i][0]>arr[i][2])
			return NULL;
		if(arr[i][1]>arr[i][3])
			return NULL;
	}
	ss_unify(arr, map, flag, size);
	for(int i=0; i<size; i++) {
		if(map[i]==i) {
			stride[i] = arr[i][2]-arr[i][0]+1;
			start[i] = arr[i][1]*fstride + arr[i][0];
			end[i] = arr[i][3]*fstride + arr[i][2];
			//if(stride[i]<=0)
			//	stride[i]=fstride;
		}
		//printf("%d-", i);
		//for(int j=0; j<4; j++)
		//printf("%d ", arr[i][j]);
		//printf("\n");
	}
	return arr;
	
}

int ss_unify(int arr[][4], long long * map, long long * flag, int size)
{
	    int is_return = 0;
	while(is_return==0) {
		is_return = 1;
		for(int i=0; i<size; i++) {
			if(map[i]!=i)
			continue;
			for(int j=i+1; j<size; j++) {
				if(map[j]!=j)
				continue;
				if(ss_isIntersect(arr, i, j)) {
					combineRect(arr, flag, i, j);
					map[j]=i;
					is_return = 0;
				}
			}
		}
	}
	for(int i=0; i<size; i++) {
		if(map[i]==i)
			continue;
		int j=i;
		while(map[j]!=j) {
			j = map[j];
		}
		map[i] = j;
		arr[i][0] -= arr[j][0];
		arr[i][1] -= arr[j][1];
	}
	
}

void combineRect(int arr[][4], long long * flag, int i, int j)
{
	arr[i][0] = ss_min(arr[i][0], arr[j][0]);
	arr[i][3] = ss_max(arr[i][3], arr[j][3]);
	arr[i][1] = ss_min(arr[i][1], arr[j][1]);
	arr[i][2] = ss_max(arr[i][2], arr[j][2]);
	flag[i] |= flag[j];
	
}

int ss_isIntersect(int arr[][4], int i, int j)
{
	if(arr[i][0]>arr[j][2] || arr[j][0]>arr[i][2])
		return 0;
	if(arr[i][1]>arr[j][3] || arr[j][1]>arr[i][3])
		return 0;
	return 1;
	
}

int ss_max(int a, int b)
{
		return (a<b) ?b :a;
}

int ss_min(int a, int b)
{
		return (a<b) ?a :b;
}

int make8ByteAligned(int num, int e_size)
{
	
	int rem = 8 - ((num*e_size)&7);
	rem = (rem==8) ?num :num+(rem/e_size);
	//printf("Align: %d %d %d\n", num, rem, e_size);
	return rem;
}

void note_down(int * p_val, int kernel_num)
{
	FILE *filePointer;
	//static int kernel_num = 0;
	char file_name[30];
	snprintf(file_name, 30, "dpu_clock_results_%d.txt", kernel_num);
	filePointer = fopen(file_name, "a");
	if(p_val==NULL) {
		fputs("\n", filePointer);
	}
	else {
		for(int i=0; i<5; i++) {
			fprintf(filePointer, "%d ", p_val[i]);
		}
		fputs("\n", filePointer);
	}
	
	//++kernel_num;
	fclose(filePointer);
	
}

int bpnn_train_kernel(BPNN * net, float * eo, float * eh)
{
	int in, hid, out;
	float out_err, hid_err;
	int sourcesize = 1024*1024;
	char * source = (char * )calloc(sourcesize, sizeof (char));
	char * kernel_bp1 = "bpnn_layerforward_ocl";
	char * kernel_bp2 = "bpnn_adjust_weights_ocl";
	char * tempchar = "./backprop_kernel.cl";
	FILE * fp = fopen(tempchar, "rb");
	int use_gpu = 1;
	cl_int err = 0;
const char * slist[2] = {source, 0};
	cl_kernel kernel1;
	cl_kernel kernel2;
	float * input_weights_one_dim;
	float * input_weights_prev_one_dim;
	float * partial_sum;
	float sum;
	float num_blocks;
	cl_program prog = CL_SUCCESS;
size_t global_work[3] = {16, 16*num_blocks, 1};
size_t local_work[3] = {16, 16, 1};
	int m = 0;
	cl_mem input_hidden_ocl;
	cl_mem input_ocl;
	cl_mem output_hidden_ocl;
	cl_mem hidden_partial_sum;
	cl_mem hidden_delta_ocl;
	cl_mem input_prev_weights_ocl;
	int _ret_val_0;
	in=net->input_n;
	hid=net->hidden_n;
	out=net->output_n;
	/* read the kernel core source */
	fread(source+strlen(source), sourcesize, 1, fp);
	fclose(fp);
	/* compile kernel */
	err=CL_SUCCESS;
	/* show warningserrors */
	/* static char log[65536]; memset(log, 0, sizeof(log)); */
	/* cl_device_id device_id = 0; */
	/* err = clGetContextInfo(context, CL_CONTEXT_DEVICES, sizeof(device_id), &device_id, NULL); */
	/* clGetProgramBuildInfo(prog, device_id, CL_PROGRAM_BUILD_LOG, sizeof(log)-1, log, NULL); */
	/* if(err || strstr(log,"warning:") || strstr(log, "error:")) printf("<<<<\n%s\n>>>>\n", log); */
	kernel1=CL_SUCCESS;
	kernel2=CL_SUCCESS;
	CL_SUCCESS;
	num_blocks=(in/16);
	input_weights_one_dim=((float * )malloc(((in+1)*(hid+1))*sizeof (float)));
	input_weights_prev_one_dim=((float * )malloc(((in+1)*(hid+1))*sizeof (float)));
	partial_sum=((float * )malloc((num_blocks*16)*sizeof (float)));
	/* set global and local workitems */
	global_work[1]=(16*num_blocks);
	/* this preprocessing stage is temporarily added to correct the bug of wrong memcopy using two-dimensional net->inputweights */
	/* todo: fix mem allocation */
	{
		int k = 0;
		#pragma cetus private(j) 
		#pragma loop name bpnn_train_kernel#0 
		for (; k<=in; k ++ )
		{
			{
				int j = 0;
				#pragma loop name bpnn_train_kernel#0#0 
				for (; j<=hid; j ++ )
				{
					input_weights_one_dim[m]=net->input_weights[k][j];
					input_weights_prev_one_dim[m]=net->input_prev_weights[k][j];
					m ++ ;
				}
			}
		}
	}
	err=CL_SUCCESS;
	err=CL_SUCCESS;
	err=CL_SUCCESS;
	err=CL_SUCCESS;
	err=CL_SUCCESS;
	err=CL_SUCCESS;
	printf("Performing GPU computation\n");
	/* write buffers */
	err=CL_SUCCESS;
	err=CL_SUCCESS;
	{
		
		struct dpu_set_t ss_dpu_set, ss_dpu;
		char* ss_profile = "backend=simulator";DPU_ASSERT(dpu_alloc(-1, ss_profile, &ss_dpu_set));
		uint32_t ss_nr_dpus = 0, ss_nr_ranks=0; 
		dpu_get_nr_dpus(ss_dpu_set, &ss_nr_dpus); 
		printf("DPU_ALLOCATED = %d\n", ss_nr_dpus); 
		dpu_get_nr_ranks(ss_dpu_set, &ss_nr_ranks); 
		printf("DPU_RANKS = %d\n", ss_nr_ranks); 
		
		
		
		
		long N_MULTI_WGS = 1;
		long Work_Group_Id = 0;
		long p_Work_Group_Id = 0;
		long bx = 0, p_bx = 0, by = 0, p_by = 0, bz = 0, p_bz = 0;
	long txl[] = {0, local_work[0]};
		long bx_limit = (global_work[0]/local_work[0]);
	long tyl[] = {0, local_work[1]};
		long by_limit = (global_work[1]/local_work[1]);
		long Work_Group_Id_Limit = (bx_limit*by_limit);
		long Thread_Id_Limit = (global_work[0]*global_work[1]);
		long N_TASKLETS = 16;
		long PARTITION_DIM_GRID = 1;
		long PARTITION_DIM_WG = 0;
		DPU_ASSERT(dpu_load(ss_dpu_set, "bpnn_layerforward_ocl", NULL));
		
		
		
		
		DPU_FOREACH(ss_dpu_set, ss_dpu)
		{
			DPU_ASSERT(dpu_copy_to(ss_dpu, "_gs", 0,  & global_work, sizeof(global_work)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "_bs", 0,  & local_work, sizeof(local_work)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "in", 0,  & in, sizeof(in)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "hid", 0,  & hid, sizeof(hid)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "dpu_gs", 1*sizeof(long),  & N_MULTI_WGS, sizeof(N_MULTI_WGS)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "GRID_LIMIT", 0,  & by_limit, sizeof(by_limit)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "N_TASKLETS", 0,  & N_TASKLETS, sizeof(N_TASKLETS)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "PARTITION_DIM_GRID", 0,  & PARTITION_DIM_GRID, sizeof(PARTITION_DIM_GRID)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "PARTITION_DIM_WG", 0,  & PARTITION_DIM_WG, sizeof(PARTITION_DIM_WG)));
		}
		DPU_ASSERT(dpu_launch(ss_dpu_set, DPU_SYNCHRONOUS));
		
		
		{
			while (Work_Group_Id<Work_Group_Id_Limit)
			{
				DPU_FOREACH(ss_dpu_set, ss_dpu)
				{
				long bxl[] = {bx, bx};
				long gxl[] = {(txl[1]*bxl[0]), ((txl[1]*bxl[1])+(txl[1]-1))};
				long byl[] = {by, ss_min(by+(N_MULTI_WGS-1), by_limit-1)};
				long gyl[] = {(tyl[1]*byl[0]), ((tyl[1]*byl[1])+(tyl[1]-1))};
				long long _bd[] = {bx, by};
					long long dpu_mram_offset = 0;
					int ss_temp_offset = 0;
					int (*ss_SE)[4];;
					long long ss_temp_length = 0;
					int n_multi_wgs = byl[1]-byl[0]+1;
					if ( ! (Work_Group_Id<Work_Group_Id_Limit))
					{
						long INPUT = -1;
						DPU_ASSERT(dpu_copy_to(ss_dpu, "INPUT", 0,  & INPUT, sizeof(INPUT)));
						continue;
					}
				long long input_ocl_map[] = {0};
				long long f_input_ocl_stride, input_ocl_offset[1] = {-1}, input_ocl_end[1] = {-1}, input_ocl_size = 0, input_ocl_cpu_offset[1] = {-1}, input_ocl_length[1] = {0}, input_ocl_ele_size;
				long long input_ocl_stride[] = {16};
				long long input_ocl_flag[] = {1};
				long long input_hidden_ocl_map[] = {0};
				long long f_input_hidden_ocl_stride, input_hidden_ocl_offset[1] = {-1}, input_hidden_ocl_end[1] = {-1}, input_hidden_ocl_size = 0, input_hidden_ocl_cpu_offset[1] = {-1}, input_hidden_ocl_length[1] = {0}, input_hidden_ocl_ele_size;
				long long input_hidden_ocl_stride[] = {1};
				long long input_hidden_ocl_flag[] = {3};
				long long hidden_partial_sum_map[] = {0};
				long long f_hidden_partial_sum_stride, hidden_partial_sum_offset[1] = {-1}, hidden_partial_sum_end[1] = {-1}, hidden_partial_sum_size = 0, hidden_partial_sum_cpu_offset[1] = {-1}, hidden_partial_sum_length[1] = {0}, hidden_partial_sum_ele_size;
				long long hidden_partial_sum_stride[] = {hid};
				long long hidden_partial_sum_flag[] = {2};
					long long input_node_size = (sizeof (float)*16);
					long long weight_matrix_size = ((sizeof (float)*16)*16);
					if(Work_Group_Id==0) {
						printf("Kernel :0 \n");
						printf("Global : ");
						printf("[0]%ld ", global_work[0]);
						printf("[1]%ld ", global_work[1]);
						printf("\n");
						printf("Local : ");
						printf("[0]%ld ", local_work[0]);
						printf("[1]%ld ", local_work[1]);
						printf("\n");
					}
					
					dpu_copy_to(ss_dpu, "_bd", 0,  & _bd, sizeof(_bd));
					DPU_ASSERT(dpu_copy_to(ss_dpu, "input_cuda_start", 0,  & dpu_mram_offset, sizeof(dpu_mram_offset)));
					input_ocl_ele_size=sizeof (float);
					if ((txl[0]<=0)&&(0<=txl[1]))
					{
					long _txl[] = {0, 0};
						input_ocl_offset[0]=(((16*byl[0])+tyl[0])+1);
						input_ocl_end[0]=(((16*byl[1])+tyl[1])+1);
					}
					f_input_ocl_stride=16;
					ss_SE = (int(*)[4])generateRectangle(input_ocl_offset, input_ocl_end, input_ocl_stride, input_ocl_map, input_ocl_flag, f_input_ocl_stride, 1);
					
					for (int copy_i = 0; copy_i<1; copy_i ++ )
					{
						input_ocl_cpu_offset[copy_i]=input_ocl_offset[copy_i];
						if (input_ocl_offset[copy_i]<0)
						{
							continue;
						}
						if (input_ocl_map[copy_i]==copy_i)
						{
							if (ss_SE!=NULL)
							{
								ss_temp_length=(((ss_SE[copy_i][3]-ss_SE[copy_i][1])+1)*input_ocl_stride[copy_i]);
							}
							else
							{
								ss_temp_length=LLONG_MAX;
							}
							if (((input_ocl_end[copy_i]-input_ocl_offset[copy_i])+1)>ss_temp_length)
							{
								int temp_ind = 0;
								float * ss_temp = ((void *)0);
								input_ocl_offset[copy_i]=0;
								input_ocl_end[copy_i]=(ss_temp_length-1);
								input_ocl_flag[copy_i]|=4;
								ss_temp_length=make8ByteAligned(ss_temp_length, sizeof(float));
								if (input_ocl_flag[copy_i]&1)
								{
									ss_temp=malloc(ss_temp_length*sizeof(float));
									for (int ss_i = ss_SE[copy_i][1]; ss_i<=ss_SE[copy_i][3]; ss_i+=1)
									{
										for (int ss_j = ss_SE[copy_i][0]; ss_j<=ss_SE[copy_i][2]; ss_j ++ )
										{
											ss_temp[temp_ind ++ ]=net->input_units[(16*ss_i)+ss_j];
										}
									}
									dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset, ss_temp, make8ByteAligned((input_ocl_end[copy_i]-input_ocl_offset[copy_i])+1, input_ocl_ele_size)*input_ocl_ele_size);
									free(ss_temp);
								}
							}
							else
							{
								int copy_to_size = (make8ByteAligned((input_ocl_end[copy_i]-input_ocl_offset[copy_i])+1, input_ocl_ele_size)*input_ocl_ele_size);
								if (input_ocl_flag[copy_i]&1)
								{
									dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset, net->input_units+input_ocl_offset[copy_i], copy_to_size);
								}
								input_ocl_stride[copy_i]=1;
							}
							input_ocl_length[copy_i]=((input_ocl_end[copy_i]-input_ocl_offset[copy_i])+1);
							input_ocl_offset[copy_i]=input_ocl_size;
							input_ocl_size=make8ByteAligned(input_ocl_length[copy_i], input_ocl_ele_size);
						}
						else
						{
							input_ocl_stride[copy_i]=input_ocl_stride[input_ocl_map[copy_i]];
							input_ocl_offset[copy_i]=(input_ocl_offset[input_ocl_map[copy_i]]+((input_ocl_stride[copy_i]*ss_SE[copy_i][1])+ss_SE[copy_i][0]));
						}
					}
					dpu_copy_to(ss_dpu, "input_cuda_offset", 0,  & input_ocl_offset, sizeof(input_ocl_offset));
					dpu_copy_to(ss_dpu, "input_cuda_stride", 0,  & input_ocl_stride, sizeof(input_ocl_stride));
					dpu_copy_to(ss_dpu, "input_cuda_flag", 0,  & input_ocl_flag, sizeof(input_ocl_flag));
					dpu_copy_to(ss_dpu, "input_cuda_length", 0,  & input_ocl_length, sizeof(input_ocl_length));
					dpu_copy_to(ss_dpu, "input_cuda_cpu_offset", 0,  & input_ocl_cpu_offset, sizeof(input_ocl_cpu_offset));
					dpu_copy_to(ss_dpu, "input_cuda_map", 0,  & input_ocl_map, sizeof(input_ocl_map));
					dpu_mram_offset+=(input_ocl_size*input_ocl_ele_size);
					
					
					DPU_ASSERT(dpu_copy_to(ss_dpu, "input_hidden_cuda_start", 0,  & dpu_mram_offset, sizeof(dpu_mram_offset)));
					input_hidden_ocl_ele_size=sizeof (float);
					input_hidden_ocl_offset[0]=(((((((hid+1)*16)*byl[0])+((hid+1)*tyl[0]))+txl[0])+1)+(hid+1));
					input_hidden_ocl_end[0]=(((((((hid+1)*16)*byl[1])+((hid+1)*tyl[1]))+txl[1])+1)+(hid+1));
					f_input_hidden_ocl_stride=1;
					ss_SE = (int(*)[4])generateRectangle(input_hidden_ocl_offset, input_hidden_ocl_end, input_hidden_ocl_stride, input_hidden_ocl_map, input_hidden_ocl_flag, f_input_hidden_ocl_stride, 1);
					
					for (int copy_i = 0; copy_i<1; copy_i ++ )
					{
						input_hidden_ocl_cpu_offset[copy_i]=input_hidden_ocl_offset[copy_i];
						if (input_hidden_ocl_offset[copy_i]<0)
						{
							continue;
						}
						if (input_hidden_ocl_map[copy_i]==copy_i)
						{
							if (ss_SE!=NULL)
							{
								ss_temp_length=(((ss_SE[copy_i][3]-ss_SE[copy_i][1])+1)*input_hidden_ocl_stride[copy_i]);
							}
							else
							{
								ss_temp_length=LLONG_MAX;
							}
							if (((input_hidden_ocl_end[copy_i]-input_hidden_ocl_offset[copy_i])+1)>ss_temp_length)
							{
								int temp_ind = 0;
								float * ss_temp = ((void *)0);
								input_hidden_ocl_offset[copy_i]=0;
								input_hidden_ocl_end[copy_i]=(ss_temp_length-1);
								input_hidden_ocl_flag[copy_i]|=4;
								ss_temp_length=make8ByteAligned(ss_temp_length, sizeof(float));
								if (input_hidden_ocl_flag[copy_i]&1)
								{
									ss_temp=malloc(ss_temp_length*sizeof(float));
									for (int ss_i = ss_SE[copy_i][1]; ss_i<=ss_SE[copy_i][3]; ss_i+=1)
									{
										for (int ss_j = ss_SE[copy_i][0]; ss_j<=ss_SE[copy_i][2]; ss_j ++ )
										{
											ss_temp[temp_ind ++ ]=input_weights_one_dim[(1*ss_i)+ss_j];
										}
									}
									dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset, ss_temp, make8ByteAligned((input_hidden_ocl_end[copy_i]-input_hidden_ocl_offset[copy_i])+1, input_hidden_ocl_ele_size)*input_hidden_ocl_ele_size);
									free(ss_temp);
								}
							}
							else
							{
								int copy_to_size = (make8ByteAligned((input_hidden_ocl_end[copy_i]-input_hidden_ocl_offset[copy_i])+1, input_hidden_ocl_ele_size)*input_hidden_ocl_ele_size);
								if (input_hidden_ocl_flag[copy_i]&1)
								{
									dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset, input_weights_one_dim+input_hidden_ocl_offset[copy_i], copy_to_size);
								}
								input_hidden_ocl_stride[copy_i]=1;
							}
							input_hidden_ocl_length[copy_i]=((input_hidden_ocl_end[copy_i]-input_hidden_ocl_offset[copy_i])+1);
							input_hidden_ocl_offset[copy_i]=input_hidden_ocl_size;
							input_hidden_ocl_size=make8ByteAligned(input_hidden_ocl_length[copy_i], input_hidden_ocl_ele_size);
						}
						else
						{
							input_hidden_ocl_stride[copy_i]=input_hidden_ocl_stride[input_hidden_ocl_map[copy_i]];
							input_hidden_ocl_offset[copy_i]=(input_hidden_ocl_offset[input_hidden_ocl_map[copy_i]]+((input_hidden_ocl_stride[copy_i]*ss_SE[copy_i][1])+ss_SE[copy_i][0]));
						}
					}
					dpu_copy_to(ss_dpu, "input_hidden_cuda_offset", 0,  & input_hidden_ocl_offset, sizeof(input_hidden_ocl_offset));
					dpu_copy_to(ss_dpu, "input_hidden_cuda_stride", 0,  & input_hidden_ocl_stride, sizeof(input_hidden_ocl_stride));
					dpu_copy_to(ss_dpu, "input_hidden_cuda_flag", 0,  & input_hidden_ocl_flag, sizeof(input_hidden_ocl_flag));
					dpu_copy_to(ss_dpu, "input_hidden_cuda_length", 0,  & input_hidden_ocl_length, sizeof(input_hidden_ocl_length));
					dpu_copy_to(ss_dpu, "input_hidden_cuda_cpu_offset", 0,  & input_hidden_ocl_cpu_offset, sizeof(input_hidden_ocl_cpu_offset));
					dpu_copy_to(ss_dpu, "input_hidden_cuda_map", 0,  & input_hidden_ocl_map, sizeof(input_hidden_ocl_map));
					dpu_mram_offset+=(input_hidden_ocl_size*input_hidden_ocl_ele_size);
					
					
					DPU_ASSERT(dpu_copy_to(ss_dpu, "hidden_partial_sum_start", 0,  & dpu_mram_offset, sizeof(dpu_mram_offset)));
					hidden_partial_sum_ele_size=sizeof (float);
					if ((txl[0]<=0)&&(0<=txl[1]))
					{
					long _txl[] = {0, 0};
						hidden_partial_sum_offset[0]=((byl[0]*hid)+tyl[0]);
						hidden_partial_sum_end[0]=((byl[1]*hid)+tyl[1]);
					}
					f_hidden_partial_sum_stride=hid;
					ss_SE = (int(*)[4])generateRectangle(hidden_partial_sum_offset, hidden_partial_sum_end, hidden_partial_sum_stride, hidden_partial_sum_map, hidden_partial_sum_flag, f_hidden_partial_sum_stride, 1);
					
					for (int copy_i = 0; copy_i<1; copy_i ++ )
					{
						hidden_partial_sum_cpu_offset[copy_i]=hidden_partial_sum_offset[copy_i];
						if (hidden_partial_sum_offset[copy_i]<0)
						{
							continue;
						}
						if (hidden_partial_sum_map[copy_i]==copy_i)
						{
							if (ss_SE!=NULL)
							{
								ss_temp_length=(((ss_SE[copy_i][3]-ss_SE[copy_i][1])+1)*hidden_partial_sum_stride[copy_i]);
							}
							else
							{
								ss_temp_length=LLONG_MAX;
							}
							if (((hidden_partial_sum_end[copy_i]-hidden_partial_sum_offset[copy_i])+1)>ss_temp_length)
							{
								int temp_ind = 0;
								hidden_partial_sum_offset[copy_i]=0;
								hidden_partial_sum_end[copy_i]=(ss_temp_length-1);
								hidden_partial_sum_flag[copy_i]|=4;
							}
							else
							{
								hidden_partial_sum_stride[copy_i]=1;
							}
							hidden_partial_sum_length[copy_i]=((hidden_partial_sum_end[copy_i]-hidden_partial_sum_offset[copy_i])+1);
							hidden_partial_sum_offset[copy_i]=hidden_partial_sum_size;
							hidden_partial_sum_size=make8ByteAligned(hidden_partial_sum_length[copy_i], hidden_partial_sum_ele_size);
						}
						else
						{
							hidden_partial_sum_stride[copy_i]=hidden_partial_sum_stride[hidden_partial_sum_map[copy_i]];
							hidden_partial_sum_offset[copy_i]=(hidden_partial_sum_offset[hidden_partial_sum_map[copy_i]]+((hidden_partial_sum_stride[copy_i]*ss_SE[copy_i][1])+ss_SE[copy_i][0]));
						}
					}
					dpu_copy_to(ss_dpu, "hidden_partial_sum_offset", 0,  & hidden_partial_sum_offset, sizeof(hidden_partial_sum_offset));
					dpu_copy_to(ss_dpu, "hidden_partial_sum_stride", 0,  & hidden_partial_sum_stride, sizeof(hidden_partial_sum_stride));
					dpu_copy_to(ss_dpu, "hidden_partial_sum_flag", 0,  & hidden_partial_sum_flag, sizeof(hidden_partial_sum_flag));
					dpu_copy_to(ss_dpu, "hidden_partial_sum_length", 0,  & hidden_partial_sum_length, sizeof(hidden_partial_sum_length));
					dpu_copy_to(ss_dpu, "hidden_partial_sum_cpu_offset", 0,  & hidden_partial_sum_cpu_offset, sizeof(hidden_partial_sum_cpu_offset));
					dpu_copy_to(ss_dpu, "hidden_partial_sum_map", 0,  & hidden_partial_sum_map, sizeof(hidden_partial_sum_map));
					dpu_mram_offset+=(hidden_partial_sum_size*hidden_partial_sum_ele_size);
					
					
					dpu_copy_to(ss_dpu, "input_node_size", 0,  & input_node_size, sizeof(input_node_size));
					dpu_copy_to(ss_dpu, "weight_matrix_size", 0,  & weight_matrix_size, sizeof(weight_matrix_size));
					// printf("WG id:%ld bxl: %ld %ld gxl: %ld %ld\n", Work_Group_Id, bxl[0], bxl[1], gxl[0], gxl[1]);
					bx=((bx+1)%bx_limit);
					if (bx==0)
					{
						by=((by+n_multi_wgs)%by_limit);
					}
					Work_Group_Id+=n_multi_wgs;
				}
				dpu_launch(ss_dpu_set, DPU_SYNCHRONOUS);
				DPU_FOREACH(ss_dpu_set, ss_dpu)
				{
				long bxl[] = {p_bx, p_bx};
				long gxl[] = {(txl[1]*bxl[0]), ((txl[1]*bxl[1])+(txl[1]-1))};
				long byl[] = {p_by, ss_min(p_by+(N_MULTI_WGS-1), by_limit-1)};
				long gyl[] = {(tyl[1]*byl[0]), ((tyl[1]*byl[1])+(tyl[1]-1))};
				long long p_bd[] = {p_bx, p_by};
					long long dpu_mram_offset = 0;
					int ss_temp_offset = 0;
					int (*ss_SE)[4];;
					long long ss_temp_length = 0;
					int n_multi_wgs = byl[1]-byl[0]+1;
					if ( ! (p_Work_Group_Id<Work_Group_Id_Limit))
					{
						break;
					}
					{
						long long f_p_hidden_partial_sum_stride, p_hidden_partial_sum_offset[1], p_hidden_partial_sum_length[1], p_hidden_partial_sum_stride[1], p_hidden_partial_sum_cpu_offset[1], p_hidden_partial_sum_flag[1], p_hidden_partial_sum_map[1], p_hidden_partial_sum_ele_size;
						DPU_ASSERT(dpu_copy_from(ss_dpu, "hidden_partial_sum_start", 0,  & dpu_mram_offset, sizeof(dpu_mram_offset)));
						p_hidden_partial_sum_ele_size=sizeof (float);
						dpu_copy_from(ss_dpu, "hidden_partial_sum_offset", 0,  & p_hidden_partial_sum_offset, sizeof(p_hidden_partial_sum_offset));
						dpu_copy_from(ss_dpu, "hidden_partial_sum_stride", 0,  & p_hidden_partial_sum_stride, sizeof(p_hidden_partial_sum_stride));
						dpu_copy_from(ss_dpu, "hidden_partial_sum_flag", 0,  & p_hidden_partial_sum_flag, sizeof(p_hidden_partial_sum_flag));
						dpu_copy_from(ss_dpu, "hidden_partial_sum_length", 0,  & p_hidden_partial_sum_length, sizeof(p_hidden_partial_sum_length));
						dpu_copy_from(ss_dpu, "hidden_partial_sum_cpu_offset", 0,  & p_hidden_partial_sum_cpu_offset, sizeof(p_hidden_partial_sum_cpu_offset));
						dpu_copy_from(ss_dpu, "hidden_partial_sum_map", 0,  & p_hidden_partial_sum_map, sizeof(p_hidden_partial_sum_map));
						f_p_hidden_partial_sum_stride=hid;
						for (int copy_i = 0; copy_i<1; copy_i ++ )
						{
							if (p_hidden_partial_sum_offset[copy_i]<0)
							{
								continue;
							}
							ss_temp_length=p_hidden_partial_sum_length[copy_i];
							if (p_hidden_partial_sum_map[copy_i]==copy_i)
							{
								if (p_hidden_partial_sum_flag[copy_i]&4)
								{
									int temp_ind = 0;
									float * ss_temp = ((void *)0);
									ss_temp_length=make8ByteAligned(ss_temp_length, sizeof(float));
									if (p_hidden_partial_sum_flag[copy_i]&2)
									{
										int ss_i_start = (p_hidden_partial_sum_cpu_offset[copy_i]/hid);
										int ss_i_limit = (ss_i_start+(p_hidden_partial_sum_length[copy_i]/p_hidden_partial_sum_stride[copy_i]));
										int ss_j_start = (p_hidden_partial_sum_cpu_offset[copy_i]%hid);
										int ss_j_limit = (p_hidden_partial_sum_stride[copy_i]+ss_j_start);
										//printf("ss_i:%d-%d, ss_j:%d-%d\n",ss_i_start, ss_i_limit, ss_j_start, ss_j_limit);
										ss_temp=malloc(ss_temp_length*sizeof(float));
										dpu_copy_from(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(p_hidden_partial_sum_offset[copy_i]*p_hidden_partial_sum_ele_size), ss_temp, make8ByteAligned(p_hidden_partial_sum_length[copy_i], p_hidden_partial_sum_ele_size)*p_hidden_partial_sum_ele_size);
										for (int ss_i = ss_i_start; ss_i<ss_i_limit; ss_i+=1)
										{
											for (int ss_j = ss_j_start; ss_j<ss_j_limit; ss_j ++ )
											{
												partial_sum[(hid*ss_i)+ss_j]=ss_temp[temp_ind ++ ];
												if (temp_ind>=p_hidden_partial_sum_length[copy_i])
												{
													break;
												}
											}
										}
										free(ss_temp);
									}
								}
								else
								{
									int copy_from_size = (make8ByteAligned(p_hidden_partial_sum_length[copy_i], p_hidden_partial_sum_ele_size)*p_hidden_partial_sum_ele_size);
									if (p_hidden_partial_sum_flag[copy_i]&2)
									{
										float * ss_temp = ((void *)0);
										int temp_ind = 0;
										if (copy_from_size>p_hidden_partial_sum_length[copy_i])
										{
											ss_temp=malloc((copy_from_size-p_hidden_partial_sum_length[copy_i])*sizeof(float));
											for (int ss_i = (p_hidden_partial_sum_cpu_offset[copy_i]+p_hidden_partial_sum_length[copy_i]); ss_i<copy_from_size; ss_i ++ )
											{
												ss_temp[temp_ind ++ ]=partial_sum[ss_i];
											}
										}
										dpu_copy_from(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(p_hidden_partial_sum_offset[copy_i]*p_hidden_partial_sum_ele_size), partial_sum+p_hidden_partial_sum_cpu_offset[copy_i], copy_from_size);
										if (copy_from_size>p_hidden_partial_sum_length[copy_i])
										{
											temp_ind=0;
											for (int ss_i = (p_hidden_partial_sum_cpu_offset[copy_i]+p_hidden_partial_sum_length[copy_i]); ss_i<copy_from_size; ss_i ++ )
											{
												partial_sum[ss_i]=ss_temp[temp_ind ++ ];
											}
											free(ss_temp);
										}
									}
								}
							}
						}
					}
					{
						long long f_p_input_hidden_ocl_stride, p_input_hidden_ocl_offset[1], p_input_hidden_ocl_length[1], p_input_hidden_ocl_stride[1], p_input_hidden_ocl_cpu_offset[1], p_input_hidden_ocl_flag[1], p_input_hidden_ocl_map[1], p_input_hidden_ocl_ele_size;
						DPU_ASSERT(dpu_copy_from(ss_dpu, "input_hidden_cuda_start", 0,  & dpu_mram_offset, sizeof(dpu_mram_offset)));
						p_input_hidden_ocl_ele_size=sizeof (float);
						dpu_copy_from(ss_dpu, "input_hidden_cuda_offset", 0,  & p_input_hidden_ocl_offset, sizeof(p_input_hidden_ocl_offset));
						dpu_copy_from(ss_dpu, "input_hidden_cuda_stride", 0,  & p_input_hidden_ocl_stride, sizeof(p_input_hidden_ocl_stride));
						dpu_copy_from(ss_dpu, "input_hidden_cuda_flag", 0,  & p_input_hidden_ocl_flag, sizeof(p_input_hidden_ocl_flag));
						dpu_copy_from(ss_dpu, "input_hidden_cuda_length", 0,  & p_input_hidden_ocl_length, sizeof(p_input_hidden_ocl_length));
						dpu_copy_from(ss_dpu, "input_hidden_cuda_cpu_offset", 0,  & p_input_hidden_ocl_cpu_offset, sizeof(p_input_hidden_ocl_cpu_offset));
						dpu_copy_from(ss_dpu, "input_hidden_cuda_map", 0,  & p_input_hidden_ocl_map, sizeof(p_input_hidden_ocl_map));
						f_p_input_hidden_ocl_stride=1;
						for (int copy_i = 0; copy_i<1; copy_i ++ )
						{
							if (p_input_hidden_ocl_offset[copy_i]<0)
							{
								continue;
							}
							ss_temp_length=p_input_hidden_ocl_length[copy_i];
							if (p_input_hidden_ocl_map[copy_i]==copy_i)
							{
								if (p_input_hidden_ocl_flag[copy_i]&4)
								{
									int temp_ind = 0;
									float * ss_temp = ((void *)0);
									ss_temp_length=make8ByteAligned(ss_temp_length, sizeof(float));
									if (p_input_hidden_ocl_flag[copy_i]&2)
									{
										int ss_i_start = (p_input_hidden_ocl_cpu_offset[copy_i]/1);
										int ss_i_limit = (ss_i_start+(p_input_hidden_ocl_length[copy_i]/p_input_hidden_ocl_stride[copy_i]));
										int ss_j_start = (p_input_hidden_ocl_cpu_offset[copy_i]%1);
										int ss_j_limit = (p_input_hidden_ocl_stride[copy_i]+ss_j_start);
										//printf("ss_i:%d-%d, ss_j:%d-%d\n",ss_i_start, ss_i_limit, ss_j_start, ss_j_limit);
										ss_temp=malloc(ss_temp_length*sizeof(float));
										dpu_copy_from(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(p_input_hidden_ocl_offset[copy_i]*p_input_hidden_ocl_ele_size), ss_temp, make8ByteAligned(p_input_hidden_ocl_length[copy_i], p_input_hidden_ocl_ele_size)*p_input_hidden_ocl_ele_size);
										for (int ss_i = ss_i_start; ss_i<ss_i_limit; ss_i+=1)
										{
											for (int ss_j = ss_j_start; ss_j<ss_j_limit; ss_j ++ )
											{
												input_weights_one_dim[(1*ss_i)+ss_j]=ss_temp[temp_ind ++ ];
												if (temp_ind>=p_input_hidden_ocl_length[copy_i])
												{
													break;
												}
											}
										}
										free(ss_temp);
									}
								}
								else
								{
									int copy_from_size = (make8ByteAligned(p_input_hidden_ocl_length[copy_i], p_input_hidden_ocl_ele_size)*p_input_hidden_ocl_ele_size);
									if (p_input_hidden_ocl_flag[copy_i]&2)
									{
										float * ss_temp = ((void *)0);
										int temp_ind = 0;
										if (copy_from_size>p_input_hidden_ocl_length[copy_i])
										{
											ss_temp=malloc((copy_from_size-p_input_hidden_ocl_length[copy_i])*sizeof(float));
											for (int ss_i = (p_input_hidden_ocl_cpu_offset[copy_i]+p_input_hidden_ocl_length[copy_i]); ss_i<copy_from_size; ss_i ++ )
											{
												ss_temp[temp_ind ++ ]=input_weights_one_dim[ss_i];
											}
										}
										dpu_copy_from(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(p_input_hidden_ocl_offset[copy_i]*p_input_hidden_ocl_ele_size), input_weights_one_dim+p_input_hidden_ocl_cpu_offset[copy_i], copy_from_size);
										if (copy_from_size>p_input_hidden_ocl_length[copy_i])
										{
											temp_ind=0;
											for (int ss_i = (p_input_hidden_ocl_cpu_offset[copy_i]+p_input_hidden_ocl_length[copy_i]); ss_i<copy_from_size; ss_i ++ )
											{
												input_weights_one_dim[ss_i]=ss_temp[temp_ind ++ ];
											}
											free(ss_temp);
										}
									}
								}
							}
						}
					}
					p_bx=((p_bx+1)%bx_limit);
					if (p_bx==0)
					{
						p_by=((p_by+n_multi_wgs)%by_limit);
					}
					p_Work_Group_Id+=n_multi_wgs;
				}
				uint64_t dpu_id = 0;
				int64_t cycles[ss_nr_dpus];
				// int64_t barrier_count[ss_nr_dpus][N_TASKLETS];
				int64_t work_done[ss_nr_dpus][N_TASKLETS];
				//int64_t T_TD_START[ss_nr_dpus][N_TASKLETS];
				//int64_t T_TD_END[ss_nr_dpus][N_TASKLETS];
				//int64_t N_WG_ID[ss_nr_dpus][N_TASKLETS];
				int ss_used_ndpus = (Work_Group_Id_Limit+N_MULTI_WGS-1)/N_MULTI_WGS;
				int ss_ndpus = (ss_nr_dpus < ss_used_ndpus) ? ss_nr_dpus : ss_used_ndpus;
				
				for(int ss_c=0; ss_c<ss_nr_dpus; ss_c++)
					cycles[ss_c]=0;
				int ss_iter_dpu;
				DPU_FOREACH(ss_dpu_set, ss_dpu, ss_iter_dpu) {
					if(ss_iter_dpu>ss_ndpus)
					break;
					//DPU_ASSERT(dpu_log_read(ss_dpu, stdout));
					DPU_ASSERT(dpu_copy_from(ss_dpu, "cycles", 0, &cycles[dpu_id], sizeof(int64_t)));
					//DPU_ASSERT(dpu_copy_from(ss_dpu, "barrier_count", 0, &barrier_count[dpu_id], sizeof(barrier_count[0])));
					DPU_ASSERT(dpu_copy_from(ss_dpu, "WORK_DONE", 0, &work_done[dpu_id], sizeof(work_done[0])));
					//DPU_ASSERT(dpu_copy_from(ss_dpu, "T_TD_START", 0, &T_TD_START[dpu_id], sizeof(T_TD_START[0])));
					//DPU_ASSERT(dpu_copy_from(ss_dpu, "T_TD_END", 0, &T_TD_END[dpu_id], sizeof(T_TD_END[0])));
					//DPU_ASSERT(dpu_copy_from(ss_dpu, "N_WG_ID", 0, &N_WG_ID[dpu_id], sizeof(N_WG_ID[0])));
					dpu_id++;
				}
				
				uint64_t total_cycles =0;
				int32_t max_count = 0;
				
				for(int i=0; i<ss_ndpus; i++) {
					printf("%ld ", cycles[i]);
					total_cycles+=cycles[i];
					if(max_count<cycles[i])
						max_count = cycles[i];
				}
				
				TOTAL_COUNT_KF += max_count;
				printf("Work_Group_Id: %ld\n", Work_Group_Id);
				printf("kernel - Total cycles = %" PRId64 "\n", total_cycles);
				printf("kernel - Max cycles = %" PRId32 "\n", max_count);
				
				printf("----------\n");
			int trace_values[5] = {N_TASKLETS, N_MULTI_WGS, PARTITION_DIM_GRID, PARTITION_DIM_WG, max_count};
				note_down(trace_values, 0);
				
				for(int i=0; i<ss_ndpus; i++) {
					for(int j=0; j<N_TASKLETS; j++) {
						//printf("%d -> T:%ld-%ld, WG:%ld, B:%ld, WD:%ld ", j, T_TD_START[i][j], T_TD_END[i][j], N_WG_ID[i][j], barrier_count[i][j], work_done[i][j]);
						if(work_done[i][j] == 0) 
								printf("[ERROR] WORK NOT DONE : Dpu %d:\n", i);
					}
					//printf("\n");
				}
				
			}
		}
		
		
		DPU_ASSERT( dpu_free(ss_dpu_set) );
		
	}
	err=CL_SUCCESS;
	{
		int j = 1;
		#pragma cetus private(k, sum) 
		#pragma loop name bpnn_train_kernel#1 
		for (; j<=hid; j ++ )
		{
			sum=0.0;
			{
				int k = 0;
				#pragma loop name bpnn_train_kernel#1#0 
				/* #pragma cetus reduction(+: sum)  */
				for (; k<num_blocks; k ++ )
				{
					sum+=partial_sum[((k*hid)+j)-1];
				}
			}
			sum+=net->input_weights[0][j];
			net->hidden_units[j]=((float)(1.0/(1.0+exp( - sum))));
		}
	}
	bpnn_layerforward(net->hidden_units, net->output_units, net->hidden_weights, hid, out);
	bpnn_output_error(net->output_delta, net->target, net->output_units, out,  & out_err);
	bpnn_hidden_error(net->hidden_delta, hid, net->output_delta, out, net->hidden_weights, net->hidden_units,  & hid_err);
	bpnn_adjust_weights(net->output_delta, out, net->hidden_units, hid, net->hidden_weights, net->hidden_prev_weights);
	err=CL_SUCCESS;
	err=CL_SUCCESS;
	err=CL_SUCCESS;
	{
		
		struct dpu_set_t ss_dpu_set, ss_dpu;
		char* ss_profile = "backend=simulator";DPU_ASSERT(dpu_alloc(-1, ss_profile, &ss_dpu_set));
		uint32_t ss_nr_dpus = 0, ss_nr_ranks=0; 
		dpu_get_nr_dpus(ss_dpu_set, &ss_nr_dpus); 
		printf("DPU_ALLOCATED = %d\n", ss_nr_dpus); 
		dpu_get_nr_ranks(ss_dpu_set, &ss_nr_ranks); 
		printf("DPU_RANKS = %d\n", ss_nr_ranks); 
		
		
		
		
		long N_MULTI_WGS = 1;
		long Work_Group_Id = 0;
		long p_Work_Group_Id = 0;
		long bx = 0, p_bx = 0, by = 0, p_by = 0, bz = 0, p_bz = 0;
	long txl[] = {0, local_work[0]};
		long bx_limit = (global_work[0]/local_work[0]);
	long tyl[] = {0, local_work[1]};
		long by_limit = (global_work[1]/local_work[1]);
		long Work_Group_Id_Limit = (bx_limit*by_limit);
		long Thread_Id_Limit = (global_work[0]*global_work[1]);
		long N_TASKLETS = 16;
		long PARTITION_DIM_GRID = 1;
		long PARTITION_DIM_WG = 0;
		DPU_ASSERT(dpu_load(ss_dpu_set, "bpnn_adjust_weights_ocl", NULL));
		
		
		
		
		DPU_FOREACH(ss_dpu_set, ss_dpu)
		{
			DPU_ASSERT(dpu_copy_to(ss_dpu, "_gs", 0,  & global_work, sizeof(global_work)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "_bs", 0,  & local_work, sizeof(local_work)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "hid", 0,  & hid, sizeof(hid)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "in", 0,  & in, sizeof(in)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "dpu_gs", 1*sizeof(long),  & N_MULTI_WGS, sizeof(N_MULTI_WGS)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "GRID_LIMIT", 0,  & by_limit, sizeof(by_limit)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "N_TASKLETS", 0,  & N_TASKLETS, sizeof(N_TASKLETS)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "PARTITION_DIM_GRID", 0,  & PARTITION_DIM_GRID, sizeof(PARTITION_DIM_GRID)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "PARTITION_DIM_WG", 0,  & PARTITION_DIM_WG, sizeof(PARTITION_DIM_WG)));
		}
		DPU_ASSERT(dpu_launch(ss_dpu_set, DPU_SYNCHRONOUS));
		
		
		{
			while (Work_Group_Id<Work_Group_Id_Limit)
			{
				DPU_FOREACH(ss_dpu_set, ss_dpu)
				{
				long bxl[] = {bx, bx};
				long gxl[] = {(txl[1]*bxl[0]), ((txl[1]*bxl[1])+(txl[1]-1))};
				long byl[] = {by, ss_min(by+(N_MULTI_WGS-1), by_limit-1)};
				long gyl[] = {(tyl[1]*byl[0]), ((tyl[1]*byl[1])+(tyl[1]-1))};
				long long _bd[] = {bx, by};
					long long dpu_mram_offset = 0;
					int ss_temp_offset = 0;
					int (*ss_SE)[4];;
					long long ss_temp_length = 0;
					int n_multi_wgs = byl[1]-byl[0]+1;
					if ( ! (Work_Group_Id<Work_Group_Id_Limit))
					{
						long INPUT = -1;
						DPU_ASSERT(dpu_copy_to(ss_dpu, "INPUT", 0,  & INPUT, sizeof(INPUT)));
						continue;
					}
				long long hidden_delta_ocl_map[] = {0};
				long long f_hidden_delta_ocl_stride, hidden_delta_ocl_offset[1] = {-1}, hidden_delta_ocl_end[1] = {-1}, hidden_delta_ocl_size = 0, hidden_delta_ocl_cpu_offset[1] = {-1}, hidden_delta_ocl_length[1] = {0}, hidden_delta_ocl_ele_size;
				long long hidden_delta_ocl_stride[] = {1};
				long long hidden_delta_ocl_flag[] = {1};
				long long input_ocl_map[] = {0};
				long long f_input_ocl_stride, input_ocl_offset[1] = {-1}, input_ocl_end[1] = {-1}, input_ocl_size = 0, input_ocl_cpu_offset[1] = {-1}, input_ocl_length[1] = {0}, input_ocl_ele_size;
				long long input_ocl_stride[] = {16};
				long long input_ocl_flag[] = {1};
				long long input_hidden_ocl_map[] = {0, 1};
				long long f_input_hidden_ocl_stride, input_hidden_ocl_offset[2] = {-1, -1}, input_hidden_ocl_end[2] = {-1, -1}, input_hidden_ocl_size = 0, input_hidden_ocl_cpu_offset[2] = {-1, -1}, input_hidden_ocl_length[2] = {0, 0}, input_hidden_ocl_ele_size;
				long long input_hidden_ocl_stride[] = {1, 1};
				long long input_hidden_ocl_flag[] = {3, 3};
				long long input_prev_weights_ocl_map[] = {0, 1};
				long long f_input_prev_weights_ocl_stride, input_prev_weights_ocl_offset[2] = {-1, -1}, input_prev_weights_ocl_end[2] = {-1, -1}, input_prev_weights_ocl_size = 0, input_prev_weights_ocl_cpu_offset[2] = {-1, -1}, input_prev_weights_ocl_length[2] = {0, 0}, input_prev_weights_ocl_ele_size;
				long long input_prev_weights_ocl_stride[] = {1, 1};
				long long input_prev_weights_ocl_flag[] = {3, 3};
				if(Work_Group_Id==0) {
					printf("Kernel :1 \n");
					printf("Global : ");
					printf("[0]%ld ", global_work[0]);
					printf("[1]%ld ", global_work[1]);
					printf("\n");
					printf("Local : ");
					printf("[0]%ld ", local_work[0]);
					printf("[1]%ld ", local_work[1]);
					printf("\n");
				}
					
					dpu_copy_to(ss_dpu, "_bd", 0,  & _bd, sizeof(_bd));
					DPU_ASSERT(dpu_copy_to(ss_dpu, "delta_start", 0,  & dpu_mram_offset, sizeof(dpu_mram_offset)));
					hidden_delta_ocl_ele_size=sizeof (float);
					hidden_delta_ocl_offset[0]=(txl[0]+1);
					hidden_delta_ocl_end[0]=(txl[1]+1);
					f_hidden_delta_ocl_stride=1;
					ss_SE = (int(*)[4])generateRectangle(hidden_delta_ocl_offset, hidden_delta_ocl_end, hidden_delta_ocl_stride, hidden_delta_ocl_map, hidden_delta_ocl_flag, f_hidden_delta_ocl_stride, 1);
					
					for (int copy_i = 0; copy_i<1; copy_i ++ )
					{
						hidden_delta_ocl_cpu_offset[copy_i]=hidden_delta_ocl_offset[copy_i];
						if (hidden_delta_ocl_offset[copy_i]<0)
						{
							continue;
						}
						if (hidden_delta_ocl_map[copy_i]==copy_i)
						{
							if (ss_SE!=NULL)
							{
								ss_temp_length=(((ss_SE[copy_i][3]-ss_SE[copy_i][1])+1)*hidden_delta_ocl_stride[copy_i]);
							}
							else
							{
								ss_temp_length=LLONG_MAX;
							}
							if (((hidden_delta_ocl_end[copy_i]-hidden_delta_ocl_offset[copy_i])+1)>ss_temp_length)
							{
								int temp_ind = 0;
								float * ss_temp = ((void *)0);
								hidden_delta_ocl_offset[copy_i]=0;
								hidden_delta_ocl_end[copy_i]=(ss_temp_length-1);
								hidden_delta_ocl_flag[copy_i]|=4;
								ss_temp_length=make8ByteAligned(ss_temp_length, sizeof(float));
								if (hidden_delta_ocl_flag[copy_i]&1)
								{
									ss_temp=malloc(ss_temp_length*sizeof(float));
									for (int ss_i = ss_SE[copy_i][1]; ss_i<=ss_SE[copy_i][3]; ss_i+=1)
									{
										for (int ss_j = ss_SE[copy_i][0]; ss_j<=ss_SE[copy_i][2]; ss_j ++ )
										{
											ss_temp[temp_ind ++ ]=net->hidden_delta[(1*ss_i)+ss_j];
										}
									}
									dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset, ss_temp, make8ByteAligned((hidden_delta_ocl_end[copy_i]-hidden_delta_ocl_offset[copy_i])+1, hidden_delta_ocl_ele_size)*hidden_delta_ocl_ele_size);
									free(ss_temp);
								}
							}
							else
							{
								int copy_to_size = (make8ByteAligned((hidden_delta_ocl_end[copy_i]-hidden_delta_ocl_offset[copy_i])+1, hidden_delta_ocl_ele_size)*hidden_delta_ocl_ele_size);
								if (hidden_delta_ocl_flag[copy_i]&1)
								{
									dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset, net->hidden_delta+hidden_delta_ocl_offset[copy_i], copy_to_size);
								}
								hidden_delta_ocl_stride[copy_i]=1;
							}
							hidden_delta_ocl_length[copy_i]=((hidden_delta_ocl_end[copy_i]-hidden_delta_ocl_offset[copy_i])+1);
							hidden_delta_ocl_offset[copy_i]=hidden_delta_ocl_size;
							hidden_delta_ocl_size=make8ByteAligned(hidden_delta_ocl_length[copy_i], hidden_delta_ocl_ele_size);
						}
						else
						{
							hidden_delta_ocl_stride[copy_i]=hidden_delta_ocl_stride[hidden_delta_ocl_map[copy_i]];
							hidden_delta_ocl_offset[copy_i]=(hidden_delta_ocl_offset[hidden_delta_ocl_map[copy_i]]+((hidden_delta_ocl_stride[copy_i]*ss_SE[copy_i][1])+ss_SE[copy_i][0]));
						}
					}
					dpu_copy_to(ss_dpu, "delta_offset", 0,  & hidden_delta_ocl_offset, sizeof(hidden_delta_ocl_offset));
					dpu_copy_to(ss_dpu, "delta_stride", 0,  & hidden_delta_ocl_stride, sizeof(hidden_delta_ocl_stride));
					dpu_copy_to(ss_dpu, "delta_flag", 0,  & hidden_delta_ocl_flag, sizeof(hidden_delta_ocl_flag));
					dpu_copy_to(ss_dpu, "delta_length", 0,  & hidden_delta_ocl_length, sizeof(hidden_delta_ocl_length));
					dpu_copy_to(ss_dpu, "delta_cpu_offset", 0,  & hidden_delta_ocl_cpu_offset, sizeof(hidden_delta_ocl_cpu_offset));
					dpu_copy_to(ss_dpu, "delta_map", 0,  & hidden_delta_ocl_map, sizeof(hidden_delta_ocl_map));
					dpu_mram_offset+=(hidden_delta_ocl_size*hidden_delta_ocl_ele_size);
					
					
					DPU_ASSERT(dpu_copy_to(ss_dpu, "ly_start", 0,  & dpu_mram_offset, sizeof(dpu_mram_offset)));
					input_ocl_ele_size=sizeof (float);
					input_ocl_offset[0]=(((16*byl[0])+tyl[0])+1);
					input_ocl_end[0]=(((16*byl[1])+tyl[1])+1);
					f_input_ocl_stride=16;
					ss_SE = (int(*)[4])generateRectangle(input_ocl_offset, input_ocl_end, input_ocl_stride, input_ocl_map, input_ocl_flag, f_input_ocl_stride, 1);
					
					for (int copy_i = 0; copy_i<1; copy_i ++ )
					{
						input_ocl_cpu_offset[copy_i]=input_ocl_offset[copy_i];
						if (input_ocl_offset[copy_i]<0)
						{
							continue;
						}
						if (input_ocl_map[copy_i]==copy_i)
						{
							if (ss_SE!=NULL)
							{
								ss_temp_length=(((ss_SE[copy_i][3]-ss_SE[copy_i][1])+1)*input_ocl_stride[copy_i]);
							}
							else
							{
								ss_temp_length=LLONG_MAX;
							}
							if (((input_ocl_end[copy_i]-input_ocl_offset[copy_i])+1)>ss_temp_length)
							{
								int temp_ind = 0;
								float * ss_temp = ((void *)0);
								input_ocl_offset[copy_i]=0;
								input_ocl_end[copy_i]=(ss_temp_length-1);
								input_ocl_flag[copy_i]|=4;
								ss_temp_length=make8ByteAligned(ss_temp_length, sizeof(float));
								if (input_ocl_flag[copy_i]&1)
								{
									ss_temp=malloc(ss_temp_length*sizeof(float));
									for (int ss_i = ss_SE[copy_i][1]; ss_i<=ss_SE[copy_i][3]; ss_i+=1)
									{
										for (int ss_j = ss_SE[copy_i][0]; ss_j<=ss_SE[copy_i][2]; ss_j ++ )
										{
											ss_temp[temp_ind ++ ]=net->input_units[(16*ss_i)+ss_j];
										}
									}
									dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset, ss_temp, make8ByteAligned((input_ocl_end[copy_i]-input_ocl_offset[copy_i])+1, input_ocl_ele_size)*input_ocl_ele_size);
									free(ss_temp);
								}
							}
							else
							{
								int copy_to_size = (make8ByteAligned((input_ocl_end[copy_i]-input_ocl_offset[copy_i])+1, input_ocl_ele_size)*input_ocl_ele_size);
								if (input_ocl_flag[copy_i]&1)
								{
									dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset, net->input_units+input_ocl_offset[copy_i], copy_to_size);
								}
								input_ocl_stride[copy_i]=1;
							}
							input_ocl_length[copy_i]=((input_ocl_end[copy_i]-input_ocl_offset[copy_i])+1);
							input_ocl_offset[copy_i]=input_ocl_size;
							input_ocl_size=make8ByteAligned(input_ocl_length[copy_i], input_ocl_ele_size);
						}
						else
						{
							input_ocl_stride[copy_i]=input_ocl_stride[input_ocl_map[copy_i]];
							input_ocl_offset[copy_i]=(input_ocl_offset[input_ocl_map[copy_i]]+((input_ocl_stride[copy_i]*ss_SE[copy_i][1])+ss_SE[copy_i][0]));
						}
					}
					dpu_copy_to(ss_dpu, "ly_offset", 0,  & input_ocl_offset, sizeof(input_ocl_offset));
					dpu_copy_to(ss_dpu, "ly_stride", 0,  & input_ocl_stride, sizeof(input_ocl_stride));
					dpu_copy_to(ss_dpu, "ly_flag", 0,  & input_ocl_flag, sizeof(input_ocl_flag));
					dpu_copy_to(ss_dpu, "ly_length", 0,  & input_ocl_length, sizeof(input_ocl_length));
					dpu_copy_to(ss_dpu, "ly_cpu_offset", 0,  & input_ocl_cpu_offset, sizeof(input_ocl_cpu_offset));
					dpu_copy_to(ss_dpu, "ly_map", 0,  & input_ocl_map, sizeof(input_ocl_map));
					dpu_mram_offset+=(input_ocl_size*input_ocl_ele_size);
					
					
					DPU_ASSERT(dpu_copy_to(ss_dpu, "w_start", 0,  & dpu_mram_offset, sizeof(dpu_mram_offset)));
					input_hidden_ocl_ele_size=sizeof (float);
					input_hidden_ocl_offset[0]=(((((((hid+1)*16)*byl[0])+((hid+1)*tyl[0]))+txl[0])+1)+(hid+1));
					input_hidden_ocl_end[0]=(((((((hid+1)*16)*byl[1])+((hid+1)*tyl[1]))+txl[1])+1)+(hid+1));
					if (((tyl[0]<=0)&&(0<=tyl[1]))&&((byl[0]<=0)&&(0<=byl[1])))
					{
					long _byl[] = {0, 0};
					long _tyl[] = {0, 0};
						input_hidden_ocl_offset[1]=(txl[0]+1);
						input_hidden_ocl_end[1]=(txl[1]+1);
					}
					f_input_hidden_ocl_stride=1;
					ss_SE = (int(*)[4])generateRectangle(input_hidden_ocl_offset, input_hidden_ocl_end, input_hidden_ocl_stride, input_hidden_ocl_map, input_hidden_ocl_flag, f_input_hidden_ocl_stride, 2);
					
					for (int copy_i = 0; copy_i<2; copy_i ++ )
					{
						input_hidden_ocl_cpu_offset[copy_i]=input_hidden_ocl_offset[copy_i];
						if (input_hidden_ocl_offset[copy_i]<0)
						{
							continue;
						}
						if (input_hidden_ocl_map[copy_i]==copy_i)
						{
							if (ss_SE!=NULL)
							{
								ss_temp_length=(((ss_SE[copy_i][3]-ss_SE[copy_i][1])+1)*input_hidden_ocl_stride[copy_i]);
							}
							else
							{
								ss_temp_length=LLONG_MAX;
							}
							if (((input_hidden_ocl_end[copy_i]-input_hidden_ocl_offset[copy_i])+1)>ss_temp_length)
							{
								int temp_ind = 0;
								float * ss_temp = ((void *)0);
								input_hidden_ocl_offset[copy_i]=0;
								input_hidden_ocl_end[copy_i]=(ss_temp_length-1);
								input_hidden_ocl_flag[copy_i]|=4;
								ss_temp_length=make8ByteAligned(ss_temp_length, sizeof(float));
								if (input_hidden_ocl_flag[copy_i]&1)
								{
									ss_temp=malloc(ss_temp_length*sizeof(float));
									for (int ss_i = ss_SE[copy_i][1]; ss_i<=ss_SE[copy_i][3]; ss_i+=1)
									{
										for (int ss_j = ss_SE[copy_i][0]; ss_j<=ss_SE[copy_i][2]; ss_j ++ )
										{
											ss_temp[temp_ind ++ ]=input_weights_one_dim[(1*ss_i)+ss_j];
										}
									}
									dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(input_hidden_ocl_size*input_hidden_ocl_ele_size), ss_temp, make8ByteAligned((input_hidden_ocl_end[copy_i]-input_hidden_ocl_offset[copy_i])+1, input_hidden_ocl_ele_size)*input_hidden_ocl_ele_size);
									free(ss_temp);
								}
							}
							else
							{
								int copy_to_size = (make8ByteAligned((input_hidden_ocl_end[copy_i]-input_hidden_ocl_offset[copy_i])+1, input_hidden_ocl_ele_size)*input_hidden_ocl_ele_size);
								if (input_hidden_ocl_flag[copy_i]&1)
								{
									dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(input_hidden_ocl_size*input_hidden_ocl_ele_size), input_weights_one_dim+input_hidden_ocl_offset[copy_i], copy_to_size);
								}
								input_hidden_ocl_stride[copy_i]=1;
							}
							input_hidden_ocl_length[copy_i]=((input_hidden_ocl_end[copy_i]-input_hidden_ocl_offset[copy_i])+1);
							input_hidden_ocl_offset[copy_i]=input_hidden_ocl_size;
							input_hidden_ocl_size+=make8ByteAligned(input_hidden_ocl_length[copy_i], input_hidden_ocl_ele_size);
						}
						else
						{
							input_hidden_ocl_stride[copy_i]=input_hidden_ocl_stride[input_hidden_ocl_map[copy_i]];
							input_hidden_ocl_offset[copy_i]=(input_hidden_ocl_offset[input_hidden_ocl_map[copy_i]]+((input_hidden_ocl_stride[copy_i]*ss_SE[copy_i][1])+ss_SE[copy_i][0]));
						}
					}
					dpu_copy_to(ss_dpu, "w_offset", 0,  & input_hidden_ocl_offset, sizeof(input_hidden_ocl_offset));
					dpu_copy_to(ss_dpu, "w_stride", 0,  & input_hidden_ocl_stride, sizeof(input_hidden_ocl_stride));
					dpu_copy_to(ss_dpu, "w_flag", 0,  & input_hidden_ocl_flag, sizeof(input_hidden_ocl_flag));
					dpu_copy_to(ss_dpu, "w_length", 0,  & input_hidden_ocl_length, sizeof(input_hidden_ocl_length));
					dpu_copy_to(ss_dpu, "w_cpu_offset", 0,  & input_hidden_ocl_cpu_offset, sizeof(input_hidden_ocl_cpu_offset));
					dpu_copy_to(ss_dpu, "w_map", 0,  & input_hidden_ocl_map, sizeof(input_hidden_ocl_map));
					dpu_mram_offset+=(input_hidden_ocl_size*input_hidden_ocl_ele_size);
					
					
					DPU_ASSERT(dpu_copy_to(ss_dpu, "oldw_start", 0,  & dpu_mram_offset, sizeof(dpu_mram_offset)));
					input_prev_weights_ocl_ele_size=sizeof (float);
					input_prev_weights_ocl_offset[0]=(((((((hid+1)*16)*byl[0])+((hid+1)*tyl[0]))+txl[0])+1)+(hid+1));
					input_prev_weights_ocl_end[0]=(((((((hid+1)*16)*byl[1])+((hid+1)*tyl[1]))+txl[1])+1)+(hid+1));
					if (((tyl[0]<=0)&&(0<=tyl[1]))&&((byl[0]<=0)&&(0<=byl[1])))
					{
					long _byl[] = {0, 0};
					long _tyl[] = {0, 0};
						input_prev_weights_ocl_offset[1]=(txl[0]+1);
						input_prev_weights_ocl_end[1]=(txl[1]+1);
					}
					f_input_prev_weights_ocl_stride=1;
					ss_SE = (int(*)[4])generateRectangle(input_prev_weights_ocl_offset, input_prev_weights_ocl_end, input_prev_weights_ocl_stride, input_prev_weights_ocl_map, input_prev_weights_ocl_flag, f_input_prev_weights_ocl_stride, 2);
					
					for (int copy_i = 0; copy_i<2; copy_i ++ )
					{
						input_prev_weights_ocl_cpu_offset[copy_i]=input_prev_weights_ocl_offset[copy_i];
						if (input_prev_weights_ocl_offset[copy_i]<0)
						{
							continue;
						}
						if (input_prev_weights_ocl_map[copy_i]==copy_i)
						{
							if (ss_SE!=NULL)
							{
								ss_temp_length=(((ss_SE[copy_i][3]-ss_SE[copy_i][1])+1)*input_prev_weights_ocl_stride[copy_i]);
							}
							else
							{
								ss_temp_length=LLONG_MAX;
							}
							if (((input_prev_weights_ocl_end[copy_i]-input_prev_weights_ocl_offset[copy_i])+1)>ss_temp_length)
							{
								int temp_ind = 0;
								float * ss_temp = ((void *)0);
								input_prev_weights_ocl_offset[copy_i]=0;
								input_prev_weights_ocl_end[copy_i]=(ss_temp_length-1);
								input_prev_weights_ocl_flag[copy_i]|=4;
								ss_temp_length=make8ByteAligned(ss_temp_length, sizeof(float));
								if (input_prev_weights_ocl_flag[copy_i]&1)
								{
									ss_temp=malloc(ss_temp_length*sizeof(float));
									for (int ss_i = ss_SE[copy_i][1]; ss_i<=ss_SE[copy_i][3]; ss_i+=1)
									{
										for (int ss_j = ss_SE[copy_i][0]; ss_j<=ss_SE[copy_i][2]; ss_j ++ )
										{
											ss_temp[temp_ind ++ ]=input_weights_prev_one_dim[(1*ss_i)+ss_j];
										}
									}
									dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(input_prev_weights_ocl_size*input_prev_weights_ocl_ele_size), ss_temp, make8ByteAligned((input_prev_weights_ocl_end[copy_i]-input_prev_weights_ocl_offset[copy_i])+1, input_prev_weights_ocl_ele_size)*input_prev_weights_ocl_ele_size);
									free(ss_temp);
								}
							}
							else
							{
								int copy_to_size = (make8ByteAligned((input_prev_weights_ocl_end[copy_i]-input_prev_weights_ocl_offset[copy_i])+1, input_prev_weights_ocl_ele_size)*input_prev_weights_ocl_ele_size);
								if (input_prev_weights_ocl_flag[copy_i]&1)
								{
									dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(input_prev_weights_ocl_size*input_prev_weights_ocl_ele_size), input_weights_prev_one_dim+input_prev_weights_ocl_offset[copy_i], copy_to_size);
								}
								input_prev_weights_ocl_stride[copy_i]=1;
							}
							input_prev_weights_ocl_length[copy_i]=((input_prev_weights_ocl_end[copy_i]-input_prev_weights_ocl_offset[copy_i])+1);
							input_prev_weights_ocl_offset[copy_i]=input_prev_weights_ocl_size;
							input_prev_weights_ocl_size+=make8ByteAligned(input_prev_weights_ocl_length[copy_i], input_prev_weights_ocl_ele_size);
						}
						else
						{
							input_prev_weights_ocl_stride[copy_i]=input_prev_weights_ocl_stride[input_prev_weights_ocl_map[copy_i]];
							input_prev_weights_ocl_offset[copy_i]=(input_prev_weights_ocl_offset[input_prev_weights_ocl_map[copy_i]]+((input_prev_weights_ocl_stride[copy_i]*ss_SE[copy_i][1])+ss_SE[copy_i][0]));
						}
					}
					dpu_copy_to(ss_dpu, "oldw_offset", 0,  & input_prev_weights_ocl_offset, sizeof(input_prev_weights_ocl_offset));
					dpu_copy_to(ss_dpu, "oldw_stride", 0,  & input_prev_weights_ocl_stride, sizeof(input_prev_weights_ocl_stride));
					dpu_copy_to(ss_dpu, "oldw_flag", 0,  & input_prev_weights_ocl_flag, sizeof(input_prev_weights_ocl_flag));
					dpu_copy_to(ss_dpu, "oldw_length", 0,  & input_prev_weights_ocl_length, sizeof(input_prev_weights_ocl_length));
					dpu_copy_to(ss_dpu, "oldw_cpu_offset", 0,  & input_prev_weights_ocl_cpu_offset, sizeof(input_prev_weights_ocl_cpu_offset));
					dpu_copy_to(ss_dpu, "oldw_map", 0,  & input_prev_weights_ocl_map, sizeof(input_prev_weights_ocl_map));
					dpu_mram_offset+=(input_prev_weights_ocl_size*input_prev_weights_ocl_ele_size);
					
					
					// printf("WG id:%ld bxl: %ld %ld gxl: %ld %ld\n", Work_Group_Id, bxl[0], bxl[1], gxl[0], gxl[1]);
					bx=((bx+1)%bx_limit);
					if (bx==0)
					{
						by=((by+n_multi_wgs)%by_limit);
					}
					Work_Group_Id+=n_multi_wgs;
				}
				dpu_launch(ss_dpu_set, DPU_SYNCHRONOUS);
				DPU_FOREACH(ss_dpu_set, ss_dpu)
				{
				long bxl[] = {p_bx, p_bx};
				long gxl[] = {(txl[1]*bxl[0]), ((txl[1]*bxl[1])+(txl[1]-1))};
				long byl[] = {p_by, ss_min(p_by+(N_MULTI_WGS-1), by_limit-1)};
				long gyl[] = {(tyl[1]*byl[0]), ((tyl[1]*byl[1])+(tyl[1]-1))};
				long long p_bd[] = {p_bx, p_by};
					long long dpu_mram_offset = 0;
					int ss_temp_offset = 0;
					int (*ss_SE)[4];;
					long long ss_temp_length = 0;
					int n_multi_wgs = byl[1]-byl[0]+1;
					if ( ! (p_Work_Group_Id<Work_Group_Id_Limit))
					{
						break;
					}
					{
						long long f_p_input_hidden_ocl_stride, p_input_hidden_ocl_offset[2], p_input_hidden_ocl_length[2], p_input_hidden_ocl_stride[2], p_input_hidden_ocl_cpu_offset[2], p_input_hidden_ocl_flag[2], p_input_hidden_ocl_map[2], p_input_hidden_ocl_ele_size;
						DPU_ASSERT(dpu_copy_from(ss_dpu, "w_start", 0,  & dpu_mram_offset, sizeof(dpu_mram_offset)));
						p_input_hidden_ocl_ele_size=sizeof (float);
						dpu_copy_from(ss_dpu, "w_offset", 0,  & p_input_hidden_ocl_offset, sizeof(p_input_hidden_ocl_offset));
						dpu_copy_from(ss_dpu, "w_stride", 0,  & p_input_hidden_ocl_stride, sizeof(p_input_hidden_ocl_stride));
						dpu_copy_from(ss_dpu, "w_flag", 0,  & p_input_hidden_ocl_flag, sizeof(p_input_hidden_ocl_flag));
						dpu_copy_from(ss_dpu, "w_length", 0,  & p_input_hidden_ocl_length, sizeof(p_input_hidden_ocl_length));
						dpu_copy_from(ss_dpu, "w_cpu_offset", 0,  & p_input_hidden_ocl_cpu_offset, sizeof(p_input_hidden_ocl_cpu_offset));
						dpu_copy_from(ss_dpu, "w_map", 0,  & p_input_hidden_ocl_map, sizeof(p_input_hidden_ocl_map));
						f_p_input_hidden_ocl_stride=1;
						for (int copy_i = 0; copy_i<2; copy_i ++ )
						{
							if (p_input_hidden_ocl_offset[copy_i]<0)
							{
								continue;
							}
							ss_temp_length=p_input_hidden_ocl_length[copy_i];
							if (p_input_hidden_ocl_map[copy_i]==copy_i)
							{
								if (p_input_hidden_ocl_flag[copy_i]&4)
								{
									int temp_ind = 0;
									float * ss_temp = ((void *)0);
									ss_temp_length=make8ByteAligned(ss_temp_length, sizeof(float));
									if (p_input_hidden_ocl_flag[copy_i]&2)
									{
										int ss_i_start = (p_input_hidden_ocl_cpu_offset[copy_i]/1);
										int ss_i_limit = (ss_i_start+(p_input_hidden_ocl_length[copy_i]/p_input_hidden_ocl_stride[copy_i]));
										int ss_j_start = (p_input_hidden_ocl_cpu_offset[copy_i]%1);
										int ss_j_limit = (p_input_hidden_ocl_stride[copy_i]+ss_j_start);
										//printf("ss_i:%d-%d, ss_j:%d-%d\n",ss_i_start, ss_i_limit, ss_j_start, ss_j_limit);
										ss_temp=malloc(ss_temp_length*sizeof(float));
										dpu_copy_from(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(p_input_hidden_ocl_offset[copy_i]*p_input_hidden_ocl_ele_size), ss_temp, make8ByteAligned(p_input_hidden_ocl_length[copy_i], p_input_hidden_ocl_ele_size)*p_input_hidden_ocl_ele_size);
										for (int ss_i = ss_i_start; ss_i<ss_i_limit; ss_i+=1)
										{
											for (int ss_j = ss_j_start; ss_j<ss_j_limit; ss_j ++ )
											{
												input_weights_one_dim[(1*ss_i)+ss_j]=ss_temp[temp_ind ++ ];
												if (temp_ind>=p_input_hidden_ocl_length[copy_i])
												{
													break;
												}
											}
										}
										free(ss_temp);
									}
								}
								else
								{
									int copy_from_size = (make8ByteAligned(p_input_hidden_ocl_length[copy_i], p_input_hidden_ocl_ele_size)*p_input_hidden_ocl_ele_size);
									if (p_input_hidden_ocl_flag[copy_i]&2)
									{
										float * ss_temp = ((void *)0);
										int temp_ind = 0;
										if (copy_from_size>p_input_hidden_ocl_length[copy_i])
										{
											ss_temp=malloc((copy_from_size-p_input_hidden_ocl_length[copy_i])*sizeof(float));
											for (int ss_i = (p_input_hidden_ocl_cpu_offset[copy_i]+p_input_hidden_ocl_length[copy_i]); ss_i<copy_from_size; ss_i ++ )
											{
												ss_temp[temp_ind ++ ]=input_weights_one_dim[ss_i];
											}
										}
										dpu_copy_from(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(p_input_hidden_ocl_offset[copy_i]*p_input_hidden_ocl_ele_size), input_weights_one_dim+p_input_hidden_ocl_cpu_offset[copy_i], copy_from_size);
										if (copy_from_size>p_input_hidden_ocl_length[copy_i])
										{
											temp_ind=0;
											for (int ss_i = (p_input_hidden_ocl_cpu_offset[copy_i]+p_input_hidden_ocl_length[copy_i]); ss_i<copy_from_size; ss_i ++ )
											{
												input_weights_one_dim[ss_i]=ss_temp[temp_ind ++ ];
											}
											free(ss_temp);
										}
									}
								}
							}
						}
					}
					{
						long long f_p_input_prev_weights_ocl_stride, p_input_prev_weights_ocl_offset[2], p_input_prev_weights_ocl_length[2], p_input_prev_weights_ocl_stride[2], p_input_prev_weights_ocl_cpu_offset[2], p_input_prev_weights_ocl_flag[2], p_input_prev_weights_ocl_map[2], p_input_prev_weights_ocl_ele_size;
						DPU_ASSERT(dpu_copy_from(ss_dpu, "oldw_start", 0,  & dpu_mram_offset, sizeof(dpu_mram_offset)));
						p_input_prev_weights_ocl_ele_size=sizeof (float);
						dpu_copy_from(ss_dpu, "oldw_offset", 0,  & p_input_prev_weights_ocl_offset, sizeof(p_input_prev_weights_ocl_offset));
						dpu_copy_from(ss_dpu, "oldw_stride", 0,  & p_input_prev_weights_ocl_stride, sizeof(p_input_prev_weights_ocl_stride));
						dpu_copy_from(ss_dpu, "oldw_flag", 0,  & p_input_prev_weights_ocl_flag, sizeof(p_input_prev_weights_ocl_flag));
						dpu_copy_from(ss_dpu, "oldw_length", 0,  & p_input_prev_weights_ocl_length, sizeof(p_input_prev_weights_ocl_length));
						dpu_copy_from(ss_dpu, "oldw_cpu_offset", 0,  & p_input_prev_weights_ocl_cpu_offset, sizeof(p_input_prev_weights_ocl_cpu_offset));
						dpu_copy_from(ss_dpu, "oldw_map", 0,  & p_input_prev_weights_ocl_map, sizeof(p_input_prev_weights_ocl_map));
						f_p_input_prev_weights_ocl_stride=1;
						for (int copy_i = 0; copy_i<2; copy_i ++ )
						{
							if (p_input_prev_weights_ocl_offset[copy_i]<0)
							{
								continue;
							}
							ss_temp_length=p_input_prev_weights_ocl_length[copy_i];
							if (p_input_prev_weights_ocl_map[copy_i]==copy_i)
							{
								if (p_input_prev_weights_ocl_flag[copy_i]&4)
								{
									int temp_ind = 0;
									float * ss_temp = ((void *)0);
									ss_temp_length=make8ByteAligned(ss_temp_length, sizeof(float));
									if (p_input_prev_weights_ocl_flag[copy_i]&2)
									{
										int ss_i_start = (p_input_prev_weights_ocl_cpu_offset[copy_i]/1);
										int ss_i_limit = (ss_i_start+(p_input_prev_weights_ocl_length[copy_i]/p_input_prev_weights_ocl_stride[copy_i]));
										int ss_j_start = (p_input_prev_weights_ocl_cpu_offset[copy_i]%1);
										int ss_j_limit = (p_input_prev_weights_ocl_stride[copy_i]+ss_j_start);
										//printf("ss_i:%d-%d, ss_j:%d-%d\n",ss_i_start, ss_i_limit, ss_j_start, ss_j_limit);
										ss_temp=malloc(ss_temp_length*sizeof(float));
										dpu_copy_from(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(p_input_prev_weights_ocl_offset[copy_i]*p_input_prev_weights_ocl_ele_size), ss_temp, make8ByteAligned(p_input_prev_weights_ocl_length[copy_i], p_input_prev_weights_ocl_ele_size)*p_input_prev_weights_ocl_ele_size);
										for (int ss_i = ss_i_start; ss_i<ss_i_limit; ss_i+=1)
										{
											for (int ss_j = ss_j_start; ss_j<ss_j_limit; ss_j ++ )
											{
												input_weights_prev_one_dim[(1*ss_i)+ss_j]=ss_temp[temp_ind ++ ];
												if (temp_ind>=p_input_prev_weights_ocl_length[copy_i])
												{
													break;
												}
											}
										}
										free(ss_temp);
									}
								}
								else
								{
									int copy_from_size = (make8ByteAligned(p_input_prev_weights_ocl_length[copy_i], p_input_prev_weights_ocl_ele_size)*p_input_prev_weights_ocl_ele_size);
									if (p_input_prev_weights_ocl_flag[copy_i]&2)
									{
										float * ss_temp = ((void *)0);
										int temp_ind = 0;
										if (copy_from_size>p_input_prev_weights_ocl_length[copy_i])
										{
											ss_temp=malloc((copy_from_size-p_input_prev_weights_ocl_length[copy_i])*sizeof(float));
											for (int ss_i = (p_input_prev_weights_ocl_cpu_offset[copy_i]+p_input_prev_weights_ocl_length[copy_i]); ss_i<copy_from_size; ss_i ++ )
											{
												ss_temp[temp_ind ++ ]=input_weights_prev_one_dim[ss_i];
											}
										}
										dpu_copy_from(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(p_input_prev_weights_ocl_offset[copy_i]*p_input_prev_weights_ocl_ele_size), input_weights_prev_one_dim+p_input_prev_weights_ocl_cpu_offset[copy_i], copy_from_size);
										if (copy_from_size>p_input_prev_weights_ocl_length[copy_i])
										{
											temp_ind=0;
											for (int ss_i = (p_input_prev_weights_ocl_cpu_offset[copy_i]+p_input_prev_weights_ocl_length[copy_i]); ss_i<copy_from_size; ss_i ++ )
											{
												input_weights_prev_one_dim[ss_i]=ss_temp[temp_ind ++ ];
											}
											free(ss_temp);
										}
									}
								}
							}
						}
					}
					p_bx=((p_bx+1)%bx_limit);
					if (p_bx==0)
					{
						p_by=((p_by+n_multi_wgs)%by_limit);
					}
					p_Work_Group_Id+=n_multi_wgs;
				}
				uint64_t dpu_id = 0;
				int64_t cycles[ss_nr_dpus];
				int64_t barrier_count[ss_nr_dpus][N_TASKLETS];
				int64_t work_done[ss_nr_dpus][N_TASKLETS];
				//int64_t T_TD_START[ss_nr_dpus][N_TASKLETS];
				//int64_t T_TD_END[ss_nr_dpus][N_TASKLETS];
				//int64_t N_WG_ID[ss_nr_dpus][N_TASKLETS];
				int ss_used_ndpus = (Work_Group_Id_Limit+N_MULTI_WGS-1)/N_MULTI_WGS;
				int ss_ndpus = (ss_nr_dpus < ss_used_ndpus) ? ss_nr_dpus : ss_used_ndpus;
				
				for(int ss_c=0; ss_c<ss_nr_dpus; ss_c++)
					cycles[ss_c]=0;
				int ss_iter_dpu;
				DPU_FOREACH(ss_dpu_set, ss_dpu, ss_iter_dpu) {
					if(ss_iter_dpu>ss_ndpus)
					break;
					//DPU_ASSERT(dpu_log_read(ss_dpu, stdout));
					DPU_ASSERT(dpu_copy_from(ss_dpu, "cycles", 0, &cycles[dpu_id], sizeof(int64_t)));
					//DPU_ASSERT(dpu_copy_from(ss_dpu, "barrier_count", 0, &barrier_count[dpu_id], sizeof(barrier_count[0])));
					DPU_ASSERT(dpu_copy_from(ss_dpu, "WORK_DONE", 0, &work_done[dpu_id], sizeof(work_done[0])));
					//DPU_ASSERT(dpu_copy_from(ss_dpu, "T_TD_START", 0, &T_TD_START[dpu_id], sizeof(T_TD_START[0])));
					//DPU_ASSERT(dpu_copy_from(ss_dpu, "T_TD_END", 0, &T_TD_END[dpu_id], sizeof(T_TD_END[0])));
					//DPU_ASSERT(dpu_copy_from(ss_dpu, "N_WG_ID", 0, &N_WG_ID[dpu_id], sizeof(N_WG_ID[0])));
					dpu_id++;
				}
				
				uint64_t total_cycles =0;
				int32_t max_count = 0;
				
				for(int i=0; i<ss_ndpus; i++) {
					printf("%ld ", cycles[i]);
					total_cycles+=cycles[i];
					if(max_count<cycles[i])
						max_count = cycles[i];
				}
				
				TOTAL_COUNT_KA += max_count;
				printf("Work_Group_Id: %ld\n", Work_Group_Id);
				printf("kernel - Total cycles = %" PRId64 "\n", total_cycles);
				printf("kernel - Max cycles = %" PRId32 "\n", max_count);
				
				printf("----------\n");
			int trace_values[5] = {N_TASKLETS, N_MULTI_WGS, PARTITION_DIM_GRID, PARTITION_DIM_WG, max_count};
				note_down(trace_values, 1);
				
				for(int i=0; i<ss_ndpus; i++) {
					for(int j=0; j<N_TASKLETS; j++) {
						//printf("%d -> T:%ld-%ld, WG:%ld, B:%ld, WD:%ld ", j, T_TD_START[i][j], T_TD_END[i][j], N_WG_ID[i][j], barrier_count[i][j], work_done[i][j]);
						if(work_done[i][j] == 0) 
								printf("[ERROR] WORK NOT DONE : Dpu %d:\n", i);
					}
					//printf("\n");
				}
				
			}
		}
		
		
		DPU_ASSERT( dpu_free(ss_dpu_set) );
		
	}
	err=CL_SUCCESS;
	err=CL_SUCCESS;
	/* newly added - samyuktha */
	err=CL_SUCCESS;
	CL_SUCCESS;
	CL_SUCCESS;
	CL_SUCCESS;
	CL_SUCCESS;
	CL_SUCCESS;
	free(input_weights_prev_one_dim);
	free(partial_sum);
	free(input_weights_one_dim);

	printf("KF: %lld, KA: %lld\n", TOTAL_COUNT_KF, TOTAL_COUNT_KA);
	return _ret_val_0;
}
