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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
/* #include <math.h> */
/* #include <iostream> */
#include <string.h>
#include "kmeans.h"
#include "open_cl.h"
/* #ifdef WIN */
/* 	#include <windows.h> */
/* #else */
/* 	#include <pthread.h> */
/* #include <systime.h> */
/* double gettime() { */
	/* 	struct timeval t; */
	/* 	gettimeofday(&t,NULL); */
	/* 	return t.tv_sec+t.tv_usec1e-6; */
/* } */
/* #endif */
/* #ifdef NV  */
/* 	#include <oclUtils.h> */
/* #else */
/* 	#include <CLcl.h> */
/* #endif */
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
	/* if( device_list ) delete device_list; */
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

cl_mem d_feature;
cl_mem d_feature_swap;
cl_mem d_cluster;
cl_mem d_membership;
cl_kernel kernel;
cl_kernel kernel_s;
cl_kernel kernel2;
int * membership_OCL;
int * membership_d;
float * feature_d;
float * clusters_d;
float * center_d;
int allocate(int npoints, int nfeatures, int nclusters, float * * feature, float * temp_feature_swap)
{
	int sourcesize = 1024*1024;
	char * source = (char * )calloc(sourcesize, sizeof (char));
	char * tempchar = "./kmeans.cl";
	FILE * fp = fopen(tempchar, "rb");
	int use_gpu = 1;
	cl_int err = 0;
const char * slist[2] = {source, 0};
	cl_program prog = CL_SUCCESS;
	char * kernel_kmeans_c = "kmeans_kernel_c";
	char * kernel_swap = "kmeans_swap";
size_t global_work[3] = {npoints, 1, 1};
size_t local_work_size[] = {256};
	int _ret_val_0;
	/* read the kernel core source */
	fread(source+strlen(source), sourcesize, 1, fp);
	fclose(fp);
	/* OpenCL initialization */
	/* compile kernel */
	err=CL_SUCCESS;
	/* show warningserrors */
	/* 	static char log[65536]; memset(log, 0, sizeof(log)); */
	/* 	cl_device_id device_id = 0; */
	/* 	err = clGetContextInfo(context, CL_CONTEXT_DEVICES, sizeof(device_id), &device_id, NULL); */
	/* 	clGetProgramBuildInfo(prog, device_id, CL_PROGRAM_BUILD_LOG, sizeof(log)-1, log, NULL); */
	/* 	if(err || strstr(log,"warning:") || strstr(log, "error:")) printf("<<<<\n%s\n>>>>\n", log); */
	kernel_s=CL_SUCCESS;
	kernel2=CL_SUCCESS;
	CL_SUCCESS;
	err=CL_SUCCESS;
	err=CL_SUCCESS;
	err=CL_SUCCESS;
	err=CL_SUCCESS;
	/* write buffers */
	err=CL_SUCCESS;
	/* Ke Wang adjustable local group size 2013/08/07 10:37:33 */
	/* work group size is defined by RD_WG_SIZE_0 or RD_WG_SIZE_0_0 201406/10 17:00:51 */
	if ((global_work[0]%local_work_size[0])!=0)
	{
		global_work[0]=(((global_work[0]/local_work_size[0])+1)*local_work_size[0]);
	}
	{
		
		struct dpu_set_t ss_dpu_set, ss_dpu;
		char* ss_profile = "backend=simulator";
		DPU_ASSERT(dpu_alloc(-1, ss_profile, &ss_dpu_set));
		uint32_t ss_nr_dpus = 0, ss_nr_ranks=0; 
		dpu_get_nr_dpus(ss_dpu_set, &ss_nr_dpus); 
		printf("DPU_ALLOCATED = %d\n", ss_nr_dpus); 
		dpu_get_nr_ranks(ss_dpu_set, &ss_nr_ranks); 
		printf("DPU_RANKS = %d\n", ss_nr_ranks); 
		
		
		
		
		long N_MULTI_WGS = 2;
		long Work_Group_Id = 0;
		long p_Work_Group_Id = 0;
		long bx = 0, p_bx = 0, by = 0, p_by = 0, bz = 0, p_bz = 0;
	long txl[] = {0, local_work_size[0]};
		long bx_limit = (global_work[0]/local_work_size[0]);
		long Work_Group_Id_Limit = bx_limit;
		long Thread_Id_Limit = global_work[0];
		long N_TASKLETS = 4;
		long PARTITION_DIM_GRID = 0;
		long PARTITION_DIM_WG = 0;
		DPU_ASSERT(dpu_load(ss_dpu_set, "kmeans_swap", NULL));
		
		
		
		
		DPU_FOREACH(ss_dpu_set, ss_dpu)
		{
			DPU_ASSERT(dpu_copy_to(ss_dpu, "_gs", 0,  & global_work, sizeof(global_work)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "_bs", 0,  & local_work_size, sizeof(local_work_size)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "npoints", 0,  & npoints, sizeof(npoints)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "nfeatures", 0,  & nfeatures, sizeof(nfeatures)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "dpu_gs", 0*sizeof(long),  & N_MULTI_WGS, sizeof(N_MULTI_WGS)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "GRID_LIMIT", 0,  & bx_limit, sizeof(bx_limit)));
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
				long bxl[] = {bx, ss_min(bx+(N_MULTI_WGS-1), bx_limit-1)};
				long gxl[] = {(txl[1]*bxl[0]), ((txl[1]*bxl[1])+(txl[1]-1))};
				long long _bd[] = {bx};
					long long dpu_mram_offset = 0;
					int ss_temp_offset = 0;
					int (*ss_SE)[4];;
					long long ss_temp_length = 0;
					int n_multi_wgs = bxl[1]-bxl[0]+1;
					if ( ! (Work_Group_Id<Work_Group_Id_Limit))
					{
						long INPUT = -1;
						DPU_ASSERT(dpu_copy_to(ss_dpu, "INPUT", 0,  & INPUT, sizeof(INPUT)));
						continue;
					}
				int host_i[] = {0, (-1+nfeatures)};
				long long d_feature_map[] = {0};
				long long f_d_feature_stride, d_feature_offset[1] = {-1}, d_feature_end[1] = {-1}, d_feature_size = 0, d_feature_cpu_offset[1] = {-1}, d_feature_length[1] = {0}, d_feature_ele_size = sizeof (float);
				long long d_feature_stride[] = {nfeatures};
				long long d_feature_flag[] = {1};
				long long d_feature_swap_map[] = {0};
				long long f_d_feature_swap_stride, d_feature_swap_offset[1] = {-1}, d_feature_swap_end[1] = {-1}, d_feature_swap_size = 0, d_feature_swap_cpu_offset[1] = {-1}, d_feature_swap_length[1] = {0}, d_feature_swap_ele_size = sizeof (float);
				long long d_feature_swap_stride[] = {npoints};
				long long d_feature_swap_flag[] = {2};
					printf("Kernel :0 \n");
					printf("Global : ");
					printf("[0]%ld ", global_work[0]);
					printf("\n");
					printf("Local : ");
					printf("[0]%ld ", local_work_size[0]);
					printf("\n");
					dpu_copy_to(ss_dpu, "_bd", 0,  & _bd, sizeof(_bd));
					DPU_ASSERT(dpu_copy_to(ss_dpu, "feature_start", 0,  & dpu_mram_offset, sizeof(dpu_mram_offset)));
					d_feature_ele_size=sizeof (float);
					if (host_i[0]<nfeatures)
					{
						if (gxl[0]<npoints)
						{
						long _host_i[] = {host_i[0], ss_min(host_i[1], nfeatures-1)};
						long _gxl[] = {gxl[0], ss_min(gxl[1], npoints-1)};
							d_feature_offset[0]=((_gxl[0]*nfeatures)+host_i[0]);
							d_feature_end[0]=((_gxl[1]*nfeatures)+host_i[1]);
						}
					}
					f_d_feature_stride=nfeatures;
					ss_SE = (int(*)[4])generateRectangle(d_feature_offset, d_feature_end, d_feature_stride, d_feature_map, d_feature_flag, f_d_feature_stride, 1);
					
					for (int copy_i = 0; copy_i<1; copy_i ++ )
					{
						d_feature_cpu_offset[copy_i]=d_feature_offset[copy_i];
						if (d_feature_offset[copy_i]<0)
						{
							continue;
						}
						if (d_feature_map[copy_i]==copy_i)
						{
							if (ss_SE!=NULL)
							{
								ss_temp_length=(((ss_SE[copy_i][3]-ss_SE[copy_i][1])+1)*d_feature_stride[copy_i]);
							}
							else
							{
								ss_temp_length=LLONG_MAX;
							}
							if (((d_feature_end[copy_i]-d_feature_offset[copy_i])+1)>ss_temp_length)
							{
								int temp_ind = 0;
								float * ss_temp = ((void *)0);
								d_feature_offset[copy_i]=0;
								d_feature_end[copy_i]=(ss_temp_length-1);
								d_feature_flag[copy_i]|=4;
								ss_temp_length=make8ByteAligned(ss_temp_length, sizeof(float));
								if (d_feature_flag[copy_i]&1)
								{
									ss_temp=malloc(ss_temp_length*sizeof(float));
									for (int ss_i = ss_SE[copy_i][1]; ss_i<=ss_SE[copy_i][3]; ss_i+=1)
									{
										for (int ss_j = ss_SE[copy_i][0]; ss_j<=ss_SE[copy_i][2]; ss_j ++ )
										{
											ss_temp[temp_ind ++ ]=feature[0][(nfeatures*ss_i)+ss_j];
										}
									}
									dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset, ss_temp, make8ByteAligned((d_feature_end[copy_i]-d_feature_offset[copy_i])+1, d_feature_ele_size)*d_feature_ele_size);
									free(ss_temp);
								}
							}
							else
							{
								int copy_to_size = (make8ByteAligned((d_feature_end[copy_i]-d_feature_offset[copy_i])+1, d_feature_ele_size)*d_feature_ele_size);
								if (d_feature_flag[copy_i]&1)
								{
									int temp_ind = 0;
									dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset, feature[0]+d_feature_offset[copy_i], copy_to_size);
								}
								d_feature_stride[copy_i]=1;
							}
							d_feature_length[copy_i]=((d_feature_end[copy_i]-d_feature_offset[copy_i])+1);
							d_feature_offset[copy_i]=d_feature_size;
							d_feature_size=make8ByteAligned(d_feature_length[copy_i], d_feature_ele_size);
						}
						else
						{
							d_feature_stride[copy_i]=d_feature_stride[d_feature_map[copy_i]];
							d_feature_offset[copy_i]=(d_feature_offset[d_feature_map[copy_i]]+((d_feature_stride[copy_i]*ss_SE[copy_i][1])+ss_SE[copy_i][0]));
						}
					}
					dpu_copy_to(ss_dpu, "feature_offset", 0,  & d_feature_offset, sizeof(d_feature_offset));
					dpu_copy_to(ss_dpu, "feature_stride", 0,  & d_feature_stride, sizeof(d_feature_stride));
					dpu_copy_to(ss_dpu, "feature_flag", 0,  & d_feature_flag, sizeof(d_feature_flag));
					dpu_copy_to(ss_dpu, "feature_length", 0,  & d_feature_length, sizeof(d_feature_length));
					dpu_copy_to(ss_dpu, "feature_cpu_offset", 0,  & d_feature_cpu_offset, sizeof(d_feature_cpu_offset));
					dpu_copy_to(ss_dpu, "feature_map", 0,  & d_feature_map, sizeof(d_feature_map));
					dpu_mram_offset+=(d_feature_size*d_feature_ele_size);
					
					
					DPU_ASSERT(dpu_copy_to(ss_dpu, "feature_swap_start", 0,  & dpu_mram_offset, sizeof(dpu_mram_offset)));
					d_feature_swap_ele_size=sizeof (float);
					if (host_i[0]<nfeatures)
					{
						if (gxl[0]<npoints)
						{
						long _host_i[] = {host_i[0], ss_min(host_i[1], nfeatures-1)};
						long _gxl[] = {gxl[0], ss_min(gxl[1], npoints-1)};
							d_feature_swap_offset[0]=((host_i[0]*npoints)+_gxl[0]);
							d_feature_swap_end[0]=((host_i[1]*npoints)+_gxl[1]);
						}
					}
					f_d_feature_swap_stride=npoints;
					ss_SE = (int(*)[4])generateRectangle(d_feature_swap_offset, d_feature_swap_end, d_feature_swap_stride, d_feature_swap_map, d_feature_swap_flag, f_d_feature_swap_stride, 1);
					
					for (int copy_i = 0; copy_i<1; copy_i ++ )
					{
						d_feature_swap_cpu_offset[copy_i]=d_feature_swap_offset[copy_i];
						if (d_feature_swap_offset[copy_i]<0)
						{
							continue;
						}
						if (d_feature_swap_map[copy_i]==copy_i)
						{
							if (ss_SE!=NULL)
							{
								ss_temp_length=(((ss_SE[copy_i][3]-ss_SE[copy_i][1])+1)*d_feature_swap_stride[copy_i]);
							}
							else
							{
								ss_temp_length=LLONG_MAX;
							}
							if (((d_feature_swap_end[copy_i]-d_feature_swap_offset[copy_i])+1)>ss_temp_length)
							{
								int temp_ind = 0;
								d_feature_swap_offset[copy_i]=0;
								d_feature_swap_end[copy_i]=(ss_temp_length-1);
								d_feature_swap_flag[copy_i]|=4;
							}
							else
							{
								d_feature_swap_stride[copy_i]=1;
							}
							d_feature_swap_length[copy_i]=((d_feature_swap_end[copy_i]-d_feature_swap_offset[copy_i])+1);
							d_feature_swap_offset[copy_i]=d_feature_swap_size;
							d_feature_swap_size=make8ByteAligned(d_feature_swap_length[copy_i], d_feature_swap_ele_size);
						}
						else
						{
							d_feature_swap_stride[copy_i]=d_feature_swap_stride[d_feature_swap_map[copy_i]];
							d_feature_swap_offset[copy_i]=(d_feature_swap_offset[d_feature_swap_map[copy_i]]+((d_feature_swap_stride[copy_i]*ss_SE[copy_i][1])+ss_SE[copy_i][0]));
						}
					}
					dpu_copy_to(ss_dpu, "feature_swap_offset", 0,  & d_feature_swap_offset, sizeof(d_feature_swap_offset));
					dpu_copy_to(ss_dpu, "feature_swap_stride", 0,  & d_feature_swap_stride, sizeof(d_feature_swap_stride));
					dpu_copy_to(ss_dpu, "feature_swap_flag", 0,  & d_feature_swap_flag, sizeof(d_feature_swap_flag));
					dpu_copy_to(ss_dpu, "feature_swap_length", 0,  & d_feature_swap_length, sizeof(d_feature_swap_length));
					dpu_copy_to(ss_dpu, "feature_swap_cpu_offset", 0,  & d_feature_swap_cpu_offset, sizeof(d_feature_swap_cpu_offset));
					dpu_copy_to(ss_dpu, "feature_swap_map", 0,  & d_feature_swap_map, sizeof(d_feature_swap_map));
					dpu_mram_offset+=(d_feature_swap_size*d_feature_swap_ele_size);
					
					
					printf("WG id:%ld bxl: %ld %ld gxl: %ld %ld\n", Work_Group_Id, bxl[0], bxl[1], gxl[0], gxl[1]);
					bx=((bx+n_multi_wgs)%bx_limit);
					Work_Group_Id+=n_multi_wgs;
				}
				dpu_launch(ss_dpu_set, DPU_SYNCHRONOUS);
				DPU_FOREACH(ss_dpu_set, ss_dpu)
				{
				long bxl[] = {p_bx, ss_min(p_bx+(N_MULTI_WGS-1), bx_limit-1)};
				long gxl[] = {(txl[1]*bxl[0]), ((txl[1]*bxl[1])+(txl[1]-1))};
				long long p_bd[] = {p_bx};
					long long dpu_mram_offset = 0;
					int ss_temp_offset = 0;
					int (*ss_SE)[4];;
					long long ss_temp_length = 0;
					int n_multi_wgs = bxl[1]-bxl[0]+1;
					if ( ! (p_Work_Group_Id<Work_Group_Id_Limit))
					{
						break;
					}
					{
						long long f_p_d_feature_swap_stride, p_d_feature_swap_offset[1], p_d_feature_swap_length[1], p_d_feature_swap_stride[1], p_d_feature_swap_cpu_offset[1], p_d_feature_swap_flag[1], p_d_feature_swap_map[1], p_d_feature_swap_ele_size = sizeof (float);
						DPU_ASSERT(dpu_copy_from(ss_dpu, "feature_swap_start", 0,  & dpu_mram_offset, sizeof(dpu_mram_offset)));
						p_d_feature_swap_ele_size=sizeof (float);
						dpu_copy_from(ss_dpu, "feature_swap_offset", 0,  & p_d_feature_swap_offset, sizeof(p_d_feature_swap_offset));
						dpu_copy_from(ss_dpu, "feature_swap_stride", 0,  & p_d_feature_swap_stride, sizeof(p_d_feature_swap_stride));
						dpu_copy_from(ss_dpu, "feature_swap_flag", 0,  & p_d_feature_swap_flag, sizeof(p_d_feature_swap_flag));
						dpu_copy_from(ss_dpu, "feature_swap_length", 0,  & p_d_feature_swap_length, sizeof(p_d_feature_swap_length));
						dpu_copy_from(ss_dpu, "feature_swap_cpu_offset", 0,  & p_d_feature_swap_cpu_offset, sizeof(p_d_feature_swap_cpu_offset));
						dpu_copy_from(ss_dpu, "feature_swap_map", 0,  & p_d_feature_swap_map, sizeof(p_d_feature_swap_map));
						f_p_d_feature_swap_stride=npoints;
						for (int copy_i = 0; copy_i<1; copy_i ++ )
						{
							if (p_d_feature_swap_offset[copy_i]<0)
							{
								continue;
							}
							ss_temp_length=p_d_feature_swap_length[copy_i];
							if (p_d_feature_swap_map[copy_i]==copy_i)
							{
								if (p_d_feature_swap_flag[copy_i]&4)
								{
									int temp_ind = 0;
									float * ss_temp = ((void *)0);
									ss_temp_length=make8ByteAligned(ss_temp_length, sizeof(float));
									if (p_d_feature_swap_flag[copy_i]&2)
									{
										int ss_i_start = (p_d_feature_swap_cpu_offset[copy_i]/npoints);
										int ss_i_limit = (ss_i_start+(p_d_feature_swap_length[copy_i]/p_d_feature_swap_stride[copy_i]));
										int ss_j_start = (p_d_feature_swap_cpu_offset[copy_i]%npoints);
										int ss_j_limit = (p_d_feature_swap_stride[copy_i]+ss_j_start);
										//printf("ss_i:%d-%d, ss_j:%d-%d\n",ss_i_start, ss_i_limit, ss_j_start, ss_j_limit);
										ss_temp=malloc(ss_temp_length*sizeof(float));
										dpu_copy_from(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(p_d_feature_swap_offset[copy_i]*p_d_feature_swap_ele_size), ss_temp, make8ByteAligned(p_d_feature_swap_length[copy_i], p_d_feature_swap_ele_size)*p_d_feature_swap_ele_size);
										for (int ss_i = ss_i_start; ss_i<ss_i_limit; ss_i+=1)
										{
											for (int ss_j = ss_j_start; ss_j<ss_j_limit; ss_j ++ )
											{
												temp_feature_swap[(npoints*ss_i)+ss_j]=ss_temp[temp_ind ++ ];
												if (temp_ind>=p_d_feature_swap_length[copy_i])
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
									int copy_from_size = (make8ByteAligned(p_d_feature_swap_length[copy_i], p_d_feature_swap_ele_size)*p_d_feature_swap_ele_size);
									if (p_d_feature_swap_flag[copy_i]&2)
									{
										float * ss_temp = ((void *)0);
										int ss_temp_length = (((ss_SE[copy_i][2]-ss_SE[copy_i][0])+1)*((ss_SE[copy_i][3]-ss_SE[copy_i][1])+1));
										int temp_ind = 0;
										if (copy_from_size>p_d_feature_swap_length[copy_i])
										{
											ss_temp=malloc((copy_from_size-p_d_feature_swap_length[copy_i])*sizeof(float));
											for (int ss_i = (p_d_feature_swap_cpu_offset[copy_i]+p_d_feature_swap_length[copy_i]); ss_i<copy_from_size; ss_i ++ )
											{
												ss_temp[temp_ind ++ ]=temp_feature_swap[ss_i];
											}
										}
										dpu_copy_from(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(p_d_feature_swap_offset[copy_i]*p_d_feature_swap_ele_size), temp_feature_swap+p_d_feature_swap_cpu_offset[copy_i], copy_from_size);
										if (copy_from_size>p_d_feature_swap_length[copy_i])
										{
											temp_ind=0;
											for (int ss_i = (p_d_feature_swap_cpu_offset[copy_i]+p_d_feature_swap_length[copy_i]); ss_i<copy_from_size; ss_i ++ )
											{
												temp_feature_swap[ss_i]=ss_temp[temp_ind ++ ];
											}
											free(ss_temp);
										}
									}
								}
							}
						}
					}
					p_bx=((p_bx+n_multi_wgs)%bx_limit);
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
	membership_OCL=((int * )malloc(npoints*sizeof (int)));
	return _ret_val_0;
}

void deallocateMemory()
{
	CL_SUCCESS;
	CL_SUCCESS;
	CL_SUCCESS;
	CL_SUCCESS;
	free(membership_OCL);
	return ;
}

int main(int argc, char * * argv)
{
	int _ret_val_0;
	printf("WG size of kernel_swap = %d, WG size of kernel_kmeans = %d \n", 256, 256);
	setup(argc, argv);
	shutdown();
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

/* in: [npoints][nfeatures] */
int kmeansOCL(float * * feature, float * temp_feature_swap, int nfeatures, int npoints, int nclusters, int * membership, float * * clusters, int * new_centers_len, float * * new_centers)
{
	int delta = 0;
	int i, j, k;
	cl_int err = 0;
size_t global_work[3] = {npoints, 1, 1};
	/* Ke Wang adjustable local group size 2013/08/07 10:37:33 */
size_t local_work_size[] = {256};
	/* work group size is defined by RD_WG_SIZE_1 or RD_WG_SIZE_1_0 201406/10 17:00:41 */
	int size = 0;
	int offset = 0;
	int _ret_val_0;
	if ((global_work[0]%local_work_size[0])!=0)
	{
		global_work[0]=(((global_work[0]/local_work_size[0])+1)*local_work_size[0]);
	}
	err=CL_SUCCESS;
	err=CL_SUCCESS;
	{
		
		struct dpu_set_t ss_dpu_set, ss_dpu;
		char* ss_profile = "backend=simulator";
		DPU_ASSERT(dpu_alloc(-1, ss_profile, &ss_dpu_set));
		uint32_t ss_nr_dpus = 0, ss_nr_ranks=0; 
		dpu_get_nr_dpus(ss_dpu_set, &ss_nr_dpus); 
		printf("DPU_ALLOCATED = %d\n", ss_nr_dpus); 
		dpu_get_nr_ranks(ss_dpu_set, &ss_nr_ranks); 
		printf("DPU_RANKS = %d\n", ss_nr_ranks); 
		
		
		
		
		long N_MULTI_WGS = 2;
		long Work_Group_Id = 0;
		long p_Work_Group_Id = 0;
		long bx = 0, p_bx = 0, by = 0, p_by = 0, bz = 0, p_bz = 0;
	long txl[] = {0, local_work_size[0]};
		long bx_limit = (global_work[0]/local_work_size[0]);
		long Work_Group_Id_Limit = bx_limit;
		long Thread_Id_Limit = global_work[0];
		long N_TASKLETS = 4;
		long PARTITION_DIM_GRID = 0;
		long PARTITION_DIM_WG = 0;
		DPU_ASSERT(dpu_load(ss_dpu_set, "kmeans_kernel_c", NULL));
		
		
		
		
		DPU_FOREACH(ss_dpu_set, ss_dpu)
		{
			DPU_ASSERT(dpu_copy_to(ss_dpu, "_gs", 0,  & global_work, sizeof(global_work)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "_bs", 0,  & local_work_size, sizeof(local_work_size)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "npoints", 0,  & npoints, sizeof(npoints)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "nclusters", 0,  & nclusters, sizeof(nclusters)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "nfeatures", 0,  & nfeatures, sizeof(nfeatures)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "offset", 0,  & offset, sizeof(offset)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "size", 0,  & size, sizeof(size)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "dpu_gs", 0*sizeof(long),  & N_MULTI_WGS, sizeof(N_MULTI_WGS)));
			DPU_ASSERT(dpu_copy_to(ss_dpu, "GRID_LIMIT", 0,  & bx_limit, sizeof(bx_limit)));
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
				long bxl[] = {bx, ss_min(bx+(N_MULTI_WGS-1), bx_limit-1)};
				long gxl[] = {(txl[1]*bxl[0]), ((txl[1]*bxl[1])+(txl[1]-1))};
				long long _bd[] = {bx};
					long long dpu_mram_offset = 0;
					int ss_temp_offset = 0;
					int (*ss_SE)[4];;
					long long ss_temp_length = 0;
					int n_multi_wgs = bxl[1]-bxl[0]+1;
					if ( ! (Work_Group_Id<Work_Group_Id_Limit))
					{
						long INPUT = -1;
						DPU_ASSERT(dpu_copy_to(ss_dpu, "INPUT", 0,  & INPUT, sizeof(INPUT)));
						continue;
					}
				int host_i[] = {0, (-1+nclusters)};
				int host_l[] = {0, (-1+nfeatures)};
				long long d_feature_swap_map[] = {0};
				long long f_d_feature_swap_stride, d_feature_swap_offset[1] = {-1}, d_feature_swap_end[1] = {-1}, d_feature_swap_size = 0, d_feature_swap_cpu_offset[1] = {-1}, d_feature_swap_length[1] = {0}, d_feature_swap_ele_size = sizeof (float);
				long long d_feature_swap_stride[] = {npoints};
				long long d_feature_swap_flag[] = {1};
				long long d_cluster_map[] = {0};
				long long f_d_cluster_stride, d_cluster_offset[1] = {-1}, d_cluster_end[1] = {-1}, d_cluster_size = 0, d_cluster_cpu_offset[1] = {-1}, d_cluster_length[1] = {0}, d_cluster_ele_size = sizeof (float);
				long long d_cluster_stride[] = {nfeatures};
				long long d_cluster_flag[] = {1};
				long long d_membership_map[] = {0};
				long long f_d_membership_stride, d_membership_offset[1] = {-1}, d_membership_end[1] = {-1}, d_membership_size = 0, d_membership_cpu_offset[1] = {-1}, d_membership_length[1] = {0}, d_membership_ele_size = sizeof (int);
				long long d_membership_stride[] = {1};
				long long d_membership_flag[] = {2};
					printf("Kernel :1 \n");
					printf("Global : ");
					printf("[0]%ld ", global_work[0]);
					printf("\n");
					printf("Local : ");
					printf("[0]%ld ", local_work_size[0]);
					printf("\n");
					dpu_copy_to(ss_dpu, "_bd", 0,  & _bd, sizeof(_bd));
					DPU_ASSERT(dpu_copy_to(ss_dpu, "feature_start", 0,  & dpu_mram_offset, sizeof(dpu_mram_offset)));
					d_feature_swap_ele_size=sizeof (float);
					if (host_l[0]<nfeatures)
					{
						if (host_i[0]<nclusters)
						{
							if (gxl[0]<npoints)
							{
							long _host_i[] = {host_i[0], ss_min(host_i[1], nclusters-1)};
							long _gxl[] = {gxl[0], ss_min(gxl[1], npoints-1)};
							long _host_l[] = {host_l[0], ss_min(host_l[1], nfeatures-1)};
								d_feature_swap_offset[0]=((host_l[0]*npoints)+_gxl[0]);
								d_feature_swap_end[0]=((host_l[1]*npoints)+_gxl[1]);
							}
						}
					}
					f_d_feature_swap_stride=npoints;
					ss_SE = (int(*)[4])generateRectangle(d_feature_swap_offset, d_feature_swap_end, d_feature_swap_stride, d_feature_swap_map, d_feature_swap_flag, f_d_feature_swap_stride, 1);
					
					for (int copy_i = 0; copy_i<1; copy_i ++ )
					{
						d_feature_swap_cpu_offset[copy_i]=d_feature_swap_offset[copy_i];
						if (d_feature_swap_offset[copy_i]<0)
						{
							continue;
						}
						if (d_feature_swap_map[copy_i]==copy_i)
						{
							if (ss_SE!=NULL)
							{
								ss_temp_length=(((ss_SE[copy_i][3]-ss_SE[copy_i][1])+1)*d_feature_swap_stride[copy_i]);
							}
							else
							{
								ss_temp_length=LLONG_MAX;
							}
							if (((d_feature_swap_end[copy_i]-d_feature_swap_offset[copy_i])+1)>ss_temp_length)
							{
								int temp_ind = 0;
								float * ss_temp = ((void *)0);
								d_feature_swap_offset[copy_i]=0;
								d_feature_swap_end[copy_i]=(ss_temp_length-1);
								d_feature_swap_flag[copy_i]|=4;
								ss_temp_length=make8ByteAligned(ss_temp_length, sizeof(float));
								if (d_feature_swap_flag[copy_i]&1)
								{
									ss_temp=malloc(ss_temp_length*sizeof(float));
									for (int ss_i = ss_SE[copy_i][1]; ss_i<=ss_SE[copy_i][3]; ss_i+=1)
									{
										for (int ss_j = ss_SE[copy_i][0]; ss_j<=ss_SE[copy_i][2]; ss_j ++ )
										{
											ss_temp[temp_ind ++ ]=temp_feature_swap[(npoints*ss_i)+ss_j];
										}
									}
									dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset, ss_temp, make8ByteAligned((d_feature_swap_end[copy_i]-d_feature_swap_offset[copy_i])+1, d_feature_swap_ele_size)*d_feature_swap_ele_size);
									free(ss_temp);
								}
							}
							else
							{
								int copy_to_size = (make8ByteAligned((d_feature_swap_end[copy_i]-d_feature_swap_offset[copy_i])+1, d_feature_swap_ele_size)*d_feature_swap_ele_size);
								if (d_feature_swap_flag[copy_i]&1)
								{
									int temp_ind = 0;
									dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset, temp_feature_swap+d_feature_swap_offset[copy_i], copy_to_size);
								}
								d_feature_swap_stride[copy_i]=1;
							}
							d_feature_swap_length[copy_i]=((d_feature_swap_end[copy_i]-d_feature_swap_offset[copy_i])+1);
							d_feature_swap_offset[copy_i]=d_feature_swap_size;
							d_feature_swap_size=make8ByteAligned(d_feature_swap_length[copy_i], d_feature_swap_ele_size);
						}
						else
						{
							d_feature_swap_stride[copy_i]=d_feature_swap_stride[d_feature_swap_map[copy_i]];
							d_feature_swap_offset[copy_i]=(d_feature_swap_offset[d_feature_swap_map[copy_i]]+((d_feature_swap_stride[copy_i]*ss_SE[copy_i][1])+ss_SE[copy_i][0]));
						}
					}
					dpu_copy_to(ss_dpu, "feature_offset", 0,  & d_feature_swap_offset, sizeof(d_feature_swap_offset));
					dpu_copy_to(ss_dpu, "feature_stride", 0,  & d_feature_swap_stride, sizeof(d_feature_swap_stride));
					dpu_copy_to(ss_dpu, "feature_flag", 0,  & d_feature_swap_flag, sizeof(d_feature_swap_flag));
					dpu_copy_to(ss_dpu, "feature_length", 0,  & d_feature_swap_length, sizeof(d_feature_swap_length));
					dpu_copy_to(ss_dpu, "feature_cpu_offset", 0,  & d_feature_swap_cpu_offset, sizeof(d_feature_swap_cpu_offset));
					dpu_copy_to(ss_dpu, "feature_map", 0,  & d_feature_swap_map, sizeof(d_feature_swap_map));
					dpu_mram_offset+=(d_feature_swap_size*d_feature_swap_ele_size);
					
					
					DPU_ASSERT(dpu_copy_to(ss_dpu, "clusters_start", 0,  & dpu_mram_offset, sizeof(dpu_mram_offset)));
					d_cluster_ele_size=sizeof (float);
					if (host_l[0]<nfeatures)
					{
						if (host_i[0]<nclusters)
						{
							if (gxl[0]<npoints)
							{
							long _host_i[] = {host_i[0], ss_min(host_i[1], nclusters-1)};
							long _gxl[] = {gxl[0], ss_min(gxl[1], npoints-1)};
							long _host_l[] = {host_l[0], ss_min(host_l[1], nfeatures-1)};
								d_cluster_offset[0]=((host_i[0]*nfeatures)+host_l[0]);
								d_cluster_end[0]=((host_i[1]*nfeatures)+host_l[1]);
							}
						}
					}
					f_d_cluster_stride=nfeatures;
					ss_SE = (int(*)[4])generateRectangle(d_cluster_offset, d_cluster_end, d_cluster_stride, d_cluster_map, d_cluster_flag, f_d_cluster_stride, 1);
					
					for (int copy_i = 0; copy_i<1; copy_i ++ )
					{
						d_cluster_cpu_offset[copy_i]=d_cluster_offset[copy_i];
						if (d_cluster_offset[copy_i]<0)
						{
							continue;
						}
						if (d_cluster_map[copy_i]==copy_i)
						{
							if (ss_SE!=NULL)
							{
								ss_temp_length=(((ss_SE[copy_i][3]-ss_SE[copy_i][1])+1)*d_cluster_stride[copy_i]);
							}
							else
							{
								ss_temp_length=LLONG_MAX;
							}
							if (((d_cluster_end[copy_i]-d_cluster_offset[copy_i])+1)>ss_temp_length)
							{
								int temp_ind = 0;
								float * ss_temp = ((void *)0);
								d_cluster_offset[copy_i]=0;
								d_cluster_end[copy_i]=(ss_temp_length-1);
								d_cluster_flag[copy_i]|=4;
								ss_temp_length=make8ByteAligned(ss_temp_length, sizeof(float));
								if (d_cluster_flag[copy_i]&1)
								{
									ss_temp=malloc(ss_temp_length*sizeof(float));
									for (int ss_i = ss_SE[copy_i][1]; ss_i<=ss_SE[copy_i][3]; ss_i+=1)
									{
										for (int ss_j = ss_SE[copy_i][0]; ss_j<=ss_SE[copy_i][2]; ss_j ++ )
										{
											ss_temp[temp_ind ++ ]=clusters[0][(nfeatures*ss_i)+ss_j];
										}
									}
									dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset, ss_temp, make8ByteAligned((d_cluster_end[copy_i]-d_cluster_offset[copy_i])+1, d_cluster_ele_size)*d_cluster_ele_size);
									free(ss_temp);
								}
							}
							else
							{
								int copy_to_size = (make8ByteAligned((d_cluster_end[copy_i]-d_cluster_offset[copy_i])+1, d_cluster_ele_size)*d_cluster_ele_size);
								if (d_cluster_flag[copy_i]&1)
								{
									int temp_ind = 0;
									dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset, clusters[0]+d_cluster_offset[copy_i], copy_to_size);
								}
								d_cluster_stride[copy_i]=1;
							}
							d_cluster_length[copy_i]=((d_cluster_end[copy_i]-d_cluster_offset[copy_i])+1);
							d_cluster_offset[copy_i]=d_cluster_size;
							d_cluster_size=make8ByteAligned(d_cluster_length[copy_i], d_cluster_ele_size);
						}
						else
						{
							d_cluster_stride[copy_i]=d_cluster_stride[d_cluster_map[copy_i]];
							d_cluster_offset[copy_i]=(d_cluster_offset[d_cluster_map[copy_i]]+((d_cluster_stride[copy_i]*ss_SE[copy_i][1])+ss_SE[copy_i][0]));
						}
					}
					dpu_copy_to(ss_dpu, "clusters_offset", 0,  & d_cluster_offset, sizeof(d_cluster_offset));
					dpu_copy_to(ss_dpu, "clusters_stride", 0,  & d_cluster_stride, sizeof(d_cluster_stride));
					dpu_copy_to(ss_dpu, "clusters_flag", 0,  & d_cluster_flag, sizeof(d_cluster_flag));
					dpu_copy_to(ss_dpu, "clusters_length", 0,  & d_cluster_length, sizeof(d_cluster_length));
					dpu_copy_to(ss_dpu, "clusters_cpu_offset", 0,  & d_cluster_cpu_offset, sizeof(d_cluster_cpu_offset));
					dpu_copy_to(ss_dpu, "clusters_map", 0,  & d_cluster_map, sizeof(d_cluster_map));
					dpu_mram_offset+=(d_cluster_size*d_cluster_ele_size);
					
					
					DPU_ASSERT(dpu_copy_to(ss_dpu, "membership_start", 0,  & dpu_mram_offset, sizeof(dpu_mram_offset)));
					d_membership_ele_size=sizeof (int);
					if (gxl[0]<npoints)
					{
					long _gxl[] = {gxl[0], ss_min(gxl[1], npoints-1)};
						d_membership_offset[0]=_gxl[0];
						d_membership_end[0]=_gxl[1];
					}
					f_d_membership_stride=1;
					ss_SE = (int(*)[4])generateRectangle(d_membership_offset, d_membership_end, d_membership_stride, d_membership_map, d_membership_flag, f_d_membership_stride, 1);
					
					for (int copy_i = 0; copy_i<1; copy_i ++ )
					{
						d_membership_cpu_offset[copy_i]=d_membership_offset[copy_i];
						if (d_membership_offset[copy_i]<0)
						{
							continue;
						}
						if (d_membership_map[copy_i]==copy_i)
						{
							if (ss_SE!=NULL)
							{
								ss_temp_length=(((ss_SE[copy_i][3]-ss_SE[copy_i][1])+1)*d_membership_stride[copy_i]);
							}
							else
							{
								ss_temp_length=LLONG_MAX;
							}
							if (((d_membership_end[copy_i]-d_membership_offset[copy_i])+1)>ss_temp_length)
							{
								int temp_ind = 0;
								d_membership_offset[copy_i]=0;
								d_membership_end[copy_i]=(ss_temp_length-1);
								d_membership_flag[copy_i]|=4;
							}
							else
							{
								d_membership_stride[copy_i]=1;
							}
							d_membership_length[copy_i]=((d_membership_end[copy_i]-d_membership_offset[copy_i])+1);
							d_membership_offset[copy_i]=d_membership_size;
							d_membership_size=make8ByteAligned(d_membership_length[copy_i], d_membership_ele_size);
						}
						else
						{
							d_membership_stride[copy_i]=d_membership_stride[d_membership_map[copy_i]];
							d_membership_offset[copy_i]=(d_membership_offset[d_membership_map[copy_i]]+((d_membership_stride[copy_i]*ss_SE[copy_i][1])+ss_SE[copy_i][0]));
						}
					}
					dpu_copy_to(ss_dpu, "membership_offset", 0,  & d_membership_offset, sizeof(d_membership_offset));
					dpu_copy_to(ss_dpu, "membership_stride", 0,  & d_membership_stride, sizeof(d_membership_stride));
					dpu_copy_to(ss_dpu, "membership_flag", 0,  & d_membership_flag, sizeof(d_membership_flag));
					dpu_copy_to(ss_dpu, "membership_length", 0,  & d_membership_length, sizeof(d_membership_length));
					dpu_copy_to(ss_dpu, "membership_cpu_offset", 0,  & d_membership_cpu_offset, sizeof(d_membership_cpu_offset));
					dpu_copy_to(ss_dpu, "membership_map", 0,  & d_membership_map, sizeof(d_membership_map));
					dpu_mram_offset+=(d_membership_size*d_membership_ele_size);
					
					
					printf("WG id:%ld bxl: %ld %ld gxl: %ld %ld\n", Work_Group_Id, bxl[0], bxl[1], gxl[0], gxl[1]);
					bx=((bx+n_multi_wgs)%bx_limit);
					Work_Group_Id+=n_multi_wgs;
				}
				dpu_launch(ss_dpu_set, DPU_SYNCHRONOUS);
				DPU_FOREACH(ss_dpu_set, ss_dpu)
				{
				long bxl[] = {p_bx, ss_min(p_bx+(N_MULTI_WGS-1), bx_limit-1)};
				long gxl[] = {(txl[1]*bxl[0]), ((txl[1]*bxl[1])+(txl[1]-1))};
				long long p_bd[] = {p_bx};
					long long dpu_mram_offset = 0;
					int ss_temp_offset = 0;
					int (*ss_SE)[4];;
					long long ss_temp_length = 0;
					int n_multi_wgs = bxl[1]-bxl[0]+1;
					if ( ! (p_Work_Group_Id<Work_Group_Id_Limit))
					{
						break;
					}
					{
						long long f_p_d_membership_stride, p_d_membership_offset[1], p_d_membership_length[1], p_d_membership_stride[1], p_d_membership_cpu_offset[1], p_d_membership_flag[1], p_d_membership_map[1], p_d_membership_ele_size = sizeof (int);
						DPU_ASSERT(dpu_copy_from(ss_dpu, "membership_start", 0,  & dpu_mram_offset, sizeof(dpu_mram_offset)));
						p_d_membership_ele_size=sizeof (int);
						dpu_copy_from(ss_dpu, "membership_offset", 0,  & p_d_membership_offset, sizeof(p_d_membership_offset));
						dpu_copy_from(ss_dpu, "membership_stride", 0,  & p_d_membership_stride, sizeof(p_d_membership_stride));
						dpu_copy_from(ss_dpu, "membership_flag", 0,  & p_d_membership_flag, sizeof(p_d_membership_flag));
						dpu_copy_from(ss_dpu, "membership_length", 0,  & p_d_membership_length, sizeof(p_d_membership_length));
						dpu_copy_from(ss_dpu, "membership_cpu_offset", 0,  & p_d_membership_cpu_offset, sizeof(p_d_membership_cpu_offset));
						dpu_copy_from(ss_dpu, "membership_map", 0,  & p_d_membership_map, sizeof(p_d_membership_map));
						f_p_d_membership_stride=1;
						for (int copy_i = 0; copy_i<1; copy_i ++ )
						{
							if (p_d_membership_offset[copy_i]<0)
							{
								continue;
							}
							ss_temp_length=p_d_membership_length[copy_i];
							if (p_d_membership_map[copy_i]==copy_i)
							{
								if (p_d_membership_flag[copy_i]&4)
								{
									int temp_ind = 0;
									int * ss_temp = ((void *)0);
									ss_temp_length=make8ByteAligned(ss_temp_length, sizeof(int));
									if (p_d_membership_flag[copy_i]&2)
									{
										int ss_i_start = (p_d_membership_cpu_offset[copy_i]/1);
										int ss_i_limit = (ss_i_start+(p_d_membership_length[copy_i]/p_d_membership_stride[copy_i]));
										int ss_j_start = (p_d_membership_cpu_offset[copy_i]%1);
										int ss_j_limit = (p_d_membership_stride[copy_i]+ss_j_start);
										//printf("ss_i:%d-%d, ss_j:%d-%d\n",ss_i_start, ss_i_limit, ss_j_start, ss_j_limit);
										ss_temp=malloc(ss_temp_length*sizeof(int));
										dpu_copy_from(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(p_d_membership_offset[copy_i]*p_d_membership_ele_size), ss_temp, make8ByteAligned(p_d_membership_length[copy_i], p_d_membership_ele_size)*p_d_membership_ele_size);
										for (int ss_i = ss_i_start; ss_i<ss_i_limit; ss_i+=1)
										{
											for (int ss_j = ss_j_start; ss_j<ss_j_limit; ss_j ++ )
											{
												membership_OCL[(1*ss_i)+ss_j]=ss_temp[temp_ind ++ ];
												if (temp_ind>=p_d_membership_length[copy_i])
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
									int copy_from_size = (make8ByteAligned(p_d_membership_length[copy_i], p_d_membership_ele_size)*p_d_membership_ele_size);
									if (p_d_membership_flag[copy_i]&2)
									{
										int * ss_temp = ((void *)0);
										int ss_temp_length = (((ss_SE[copy_i][2]-ss_SE[copy_i][0])+1)*((ss_SE[copy_i][3]-ss_SE[copy_i][1])+1));
										int temp_ind = 0;
										if (copy_from_size>p_d_membership_length[copy_i])
										{
											ss_temp=malloc((copy_from_size-p_d_membership_length[copy_i])*sizeof(int));
											for (int ss_i = (p_d_membership_cpu_offset[copy_i]+p_d_membership_length[copy_i]); ss_i<copy_from_size; ss_i ++ )
											{
												ss_temp[temp_ind ++ ]=membership_OCL[ss_i];
											}
										}
										dpu_copy_from(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(p_d_membership_offset[copy_i]*p_d_membership_ele_size), membership_OCL+p_d_membership_cpu_offset[copy_i], copy_from_size);
										if (copy_from_size>p_d_membership_length[copy_i])
										{
											temp_ind=0;
											for (int ss_i = (p_d_membership_cpu_offset[copy_i]+p_d_membership_length[copy_i]); ss_i<copy_from_size; ss_i ++ )
											{
												membership_OCL[ss_i]=ss_temp[temp_ind ++ ];
											}
											free(ss_temp);
										}
									}
								}
							}
						}
					}
					p_bx=((p_bx+n_multi_wgs)%bx_limit);
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
	CL_SUCCESS;
	err=CL_SUCCESS;
	delta=0;
	#pragma cetus private(cluster_id, j) 
	#pragma loop name kmeansOCL#0 
	/* #pragma cetus reduction(+: delta, new_centers[cluster_id][j])  */
	for (i=0; i<npoints; i ++ )
	{
		int cluster_id = membership_OCL[i];
		new_centers_len[cluster_id] ++ ;
		if (membership_OCL[i]!=membership[i])
		{
			delta ++ ;
			membership[i]=membership_OCL[i];
		}
		#pragma cetus private(j) 
		#pragma loop name kmeansOCL#0#0 
		for (j=0; j<nfeatures; j ++ )
		{
			new_centers[cluster_id][j]+=feature[i][j];
		}
	}
	return delta;
}
