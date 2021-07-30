#include<dpu.h>
#include<inttypes.h>
#include<math.h>
void note_down(int* p_val, int kernel_num, int t);
int make8ByteAligned(int num, int e_size);
int ss_min(int a, int b);
int ss_max(int a, int b);
int* generateRectangle(long long *start, long long *end, long long *stride, long long *map, long long *flag, int fstride, int size);
int ss_unify(int arr[][4], long long * map, long long* flag, int size);
void combineRect(int arr[][4], long long *flag, int i, int j);
int ss_isIntersect(int arr[][4], int i, int j);

int n_tasklets=0;
int group_size=0;

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
#include "gaussianElim.h"
/* #include <math.h> */

#ifdef RD_WG_SIZE_0_0

        #define BLOCK_SIZE_0 RD_WG_SIZE_0_0

#elif defined(RD_WG_SIZE_0)

        #define BLOCK_SIZE_0 RD_WG_SIZE_0

#elif defined(RD_WG_SIZE)

        #define BLOCK_SIZE_0 RD_WG_SIZE

#else

        #define BLOCK_SIZE_0 0

#endif



                                                

#ifdef RD_WG_SIZE_1_0

        #define BLOCK_SIZE_1_X RD_WG_SIZE_1_0

#elif defined(RD_WG_SIZE_1)

        #define BLOCK_SIZE_1_X RD_WG_SIZE_1

#elif defined(RD_WG_SIZE)

        #define BLOCK_SIZE_1_X RD_WG_SIZE

#else

        #define BLOCK_SIZE_1_X 0

#endif



#ifdef RD_WG_SIZE_1_1

        #define BLOCK_SIZE_1_Y RD_WG_SIZE_1_1

#elif defined(RD_WG_SIZE_1)

        #define BLOCK_SIZE_1_Y RD_WG_SIZE_1

#elif defined(RD_WG_SIZE)

        #define BLOCK_SIZE_1_Y RD_WG_SIZE

#else

        #define BLOCK_SIZE_1_Y 0

#endif

cl_context context = (void * )0;
/* create both matrix and right hand side, Ke Wang 201308/12 11:51:06 */
void create_matrix(float * m, int size)
{
	int i, j;
	float lamda =  - 0.01;
	float coe[((2*size)-1)];
	float coe_i = 0.0;
	#pragma loop name create_matrix#0 
	for (i=0; i<size; i ++ )
	{
		coe_i=(10*exp(lamda*i));
		j=((size-1)+i);
		coe[j]=coe_i;
		j=((size-1)-i);
		coe[j]=coe_i;
	}
	#pragma loop name create_matrix#1 
	for (i=0; i<size; i ++ )
	{
		#pragma loop name create_matrix#1#0 
		for (j=0; j<size; j ++ )
		{
			m[(i*size)+j]=coe[((size-1)-i)+j];
		}
	}
	return ;
}

int main(int argc, char * argv[])
{
	float * a = (void * )0, * b = (void * )0, * finalVec = (void * )0;
	float * m = (void * )0;
	int size =  - 1;
	FILE * fp;
	char filename[200];
	int quiet = 0, timing = 0, platform =  - 1, device =  - 1;
	int _ret_val_0;
	printf("WG size of kernel 1 = %d, WG size of kernel 2= %d X %d\n", BLOCK_SIZE_0, BLOCK_SIZE_1_X, BLOCK_SIZE_1_Y);
	/* args */
	/* parse command line */
	if (parseCommandline(argc, argv, filename,  & quiet,  & timing,  & platform,  & device,  & size))
	{
		printUsage();
		_ret_val_0=0;
		return _ret_val_0;
	}
	context=CL_SUCCESS;
	if (size<1)
	{
		fp=fopen(filename, "r");
		fscanf(fp, "%d",  & size);
		a=((float * )malloc((size*size)*sizeof (float)));
		InitMat(fp, size, a, size, size);
		b=((float * )malloc(size*sizeof (float)));
		InitAry(fp, b, size);
		fclose(fp);
	}
	else
	{
		printf("create input internally before create, size = %d \n", size);
		a=((float * )malloc((size*size)*sizeof (float)));
		create_matrix(a, size);
		b=((float * )malloc(size*sizeof (float)));
		{
			int i = 0;
			#pragma loop name main#0 
			for (; i<size; i ++ )
			{
				b[i]=1.0;
			}
		}
	}
	// if ( ! quiet)
	// {
	// 	printf("The input matrix a is:\n");
	// 	PrintMat(a, size, size, size);
	// 	printf("The input array b is:\n");
	// 	PrintAry(b, size);
	// }
	/* create the solution matrix */
	m=((float * )malloc((size*size)*sizeof (float)));
	/* create a new vector to hold the final answer */
	finalVec=((float * )malloc(size*sizeof (float)));
	InitPerRun(size, m);
	/* begin timing	 */
	/* printf("The result of array b is before run: \n"); */
	/* PrintAry(b, size); */
	/* run kernels */
	printf("G%d NT%d\n", group_size, n_tasklets);
	ForwardSub(context, a, b, m, size, timing);
	/* printf("The result of array b is after run: \n"); */
	/* PrintAry(b, size); */
	/* end timing */
	if ( ! quiet)
	{
		// printf("The result of matrix m is: \n");
		// PrintMat(m, size, size, size);
		// printf("The result of matrix a is: \n");
		// PrintMat(a, size, size, size);
		printf("The result of array b is: \n");
		PrintAry(b, size);
		BackSub(a, b, finalVec, size);
		printf("The final solution is: \n");
		PrintAry(finalVec, size);
	}
	free(m);
	free(a);
	free(b);
	free(finalVec);
	CL_SUCCESS;
	/* OpenClGaussianElimination(context,timing); */
	_ret_val_0=0;
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
		}
		// printf("%d-", i);
		// for(int j=0; j<4; j++)
		// printf("%d ", arr[i][j]);
		// printf("\n");
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
	// printf("Align: %d %d %d\n", num, rem, e_size);
	return rem;
}

void note_down(int * p_val, int kernel_num, int t)
{
	FILE *filePointer;
	// static int kernel_num = 0;
	char file_name[50];
	snprintf(file_name, 50, "dpu_clock_results_%d_mat208_all.txt", kernel_num);
	filePointer = fopen(file_name, "a");
	if(p_val==NULL) {
		fputs("\n", filePointer);
	}
	else {
		fprintf(filePointer, "%d ", t);
		for(int i=0; i<5; i++) {
			fprintf(filePointer, "%d ", p_val[i]);
		}

		fputs("\n", filePointer);
	}
	
	// ++kernel_num;
	fclose(filePointer);
	
}


/*
------------------------------------------------------

 ForwardSub() -- Forward substitution of Gaussian

 ** elimination.

 **------------------------------------------------------


*/
void ForwardSub(cl_context context, float * a, float * b, float * m, long long size, int timing)
{
	/* NOT ORIGINAL */
	/* 1. set up kernels */
	cl_kernel fan1_kernel, fan2_kernel;
	cl_int status = 0;
	cl_program gaussianElim_program;
	cl_event writeEvent, kernelEvent, readEvent;
	float writeTime = 0, readTime = 0, kernelTime = 0;
	float writeMB = 0, readMB = 0;
	cl_mem a_dev, b_dev, m_dev;
	cl_int error = 0;
	cl_command_queue command_queue = CL_SUCCESS;
	size_t globalWorksizeFan1[1];
	size_t globalWorksizeFan2[2];
size_t localWorksizeFan1Buf[1] = {BLOCK_SIZE_0};
size_t localWorksizeFan2Buf[2] = {BLOCK_SIZE_1_X, BLOCK_SIZE_1_Y};
	size_t * localWorksizeFan1 = (void * )0;
	size_t * localWorksizeFan2 = (void * )0;
	long t;
	gaussianElim_program=CL_SUCCESS;
	fan1_kernel=CL_SUCCESS;
	status=CL_SUCCESS;
	if (status)
	{
		exit(1);
	}
	fan2_kernel=CL_SUCCESS;
	status=CL_SUCCESS;
	if (status)
	{
		exit(1);
	}
	/* 2. set up memory on device and send ipts data to device */
	error=CL_SUCCESS;
	error=CL_SUCCESS;
	error=CL_SUCCESS;
	/* change to 0 for nonblocking write */
	/* offset */
	error=CL_SUCCESS;
	if (timing)
	{
		writeTime+=eventTime(writeEvent, command_queue);
	}
	CL_SUCCESS;
	/* change to 0 for nonblocking write */
	/* offset */
	error=CL_SUCCESS;
	if (timing)
	{
		writeTime+=eventTime(writeEvent, command_queue);
	}
	CL_SUCCESS;
	/* change to 0 for nonblocking write */
	/* offset */
	error=CL_SUCCESS;
	if (timing)
	{
		writeTime+=eventTime(writeEvent, command_queue);
	}
	CL_SUCCESS;
	writeMB=((float)(((sizeof (float)*size)*((size+size)+1))/1000000.0));
	/* 3. Determine block sizes */
	globalWorksizeFan1[0]=size;
	globalWorksizeFan2[0]=size;
	globalWorksizeFan2[1]=size;
	if (localWorksizeFan1Buf[0])
	{
		localWorksizeFan1=localWorksizeFan1Buf;
		globalWorksizeFan1[0]=(((int)ceil(globalWorksizeFan1[0]/((double)localWorksizeFan1Buf[0])))*localWorksizeFan1Buf[0]);
	}
	if (localWorksizeFan2Buf[0])
	{
		localWorksizeFan2=localWorksizeFan2Buf;
		globalWorksizeFan2[0]=(((int)ceil(globalWorksizeFan2[0]/((double)localWorksizeFan2Buf[0])))*localWorksizeFan2Buf[0]);
		globalWorksizeFan2[1]=(((int)ceil(globalWorksizeFan2[1]/((double)localWorksizeFan2Buf[1])))*localWorksizeFan2Buf[1]);
	}
	/* 4. Setup and Run kernels */
	#pragma loop name ForwardSub#0 
	/* #pragma cetus reduction(+: kernelTime)  */
	long long FINAL_TOTAL_CYCLES_0 = 0;
	long long FINAL_TOTAL_CYCLES_1 = 0;

	FILE *filePointer1 = fopen("mat_208_k1.txt", "a");
	FILE *filePointer2 = fopen("mat_208_k2.txt", "a");
	FILE *filePointer3 = fopen("mat_208_k12.txt", "a");
	


	for (t=0; t<(size-1); t ++ )
	{
		/* kernel args */
		cl_int argchk;
		argchk=CL_SUCCESS;
		argchk|=CL_SUCCESS;
		argchk|=CL_SUCCESS;
		argchk|=CL_SUCCESS;
		argchk|=CL_SUCCESS;
		CL_SUCCESS;
		/* launch kernel */
		{
			
			struct dpu_set_t ss_dpu_set, ss_dpu;
			char* ss_profile = "backend=simulator";DPU_ASSERT(dpu_alloc(-1, ss_profile, &ss_dpu_set));
			uint32_t ss_nr_dpus = 0, ss_nr_ranks=0; 
			dpu_get_nr_dpus(ss_dpu_set, &ss_nr_dpus); 

			
			
			
			
			
			long long N_MULTI_WGS = group_size;
			long Work_Group_Id = 0;
			long p_Work_Group_Id = 0;
			long bx = 0, p_bx = 0, by = 0, p_by = 0, bz = 0, p_bz = 0;
		long txl[] = {0, localWorksizeFan1[0]};
			long bx_limit = (globalWorksizeFan1[0]/localWorksizeFan1[0]);
			long Work_Group_Id_Limit = bx_limit;
			long Thread_Id_Limit = globalWorksizeFan1[0];
			long N_TASKLETS = n_tasklets;
			long PARTITION_DIM_GRID = 0;
			long PARTITION_DIM_WG = 0;
			DPU_ASSERT(dpu_load(ss_dpu_set, "Fan1", NULL));

			if(t==size-2) {
			printf("DPU_ALLOCATED = %d\n", ss_nr_dpus); 
			dpu_get_nr_ranks(ss_dpu_set, &ss_nr_ranks); 
			printf("DPU_RANKS = %d\n", ss_nr_ranks); 
			printf("Work_Group_Id_Limit = %ld\n", Work_Group_Id_Limit);
						printf("Kernel :0 \n");
						printf("Global : ");
						printf("[0]%ld ", globalWorksizeFan1[0]);
						printf("\n");
						printf("Local : ");
						printf("[0]%ld ", localWorksizeFan1[0]);
						printf("\n");

			}
			
			
			
			
			DPU_FOREACH(ss_dpu_set, ss_dpu)
			{
				DPU_ASSERT(dpu_copy_to(ss_dpu, "GRID_LIMIT", 0, &bx_limit, sizeof(bx_limit)));
				DPU_ASSERT(dpu_copy_to(ss_dpu, "_gs", 0, globalWorksizeFan1, sizeof(globalWorksizeFan1)));
				DPU_ASSERT(dpu_copy_to(ss_dpu, "_bs", 0, localWorksizeFan1, sizeof(localWorksizeFan1)));
				DPU_ASSERT(dpu_copy_to(ss_dpu, "size", 0,  & size, sizeof(size)));
				DPU_ASSERT(dpu_copy_to(ss_dpu, "t", 0,  & t, sizeof(t)));
				DPU_ASSERT(dpu_copy_to(ss_dpu, "dpu_gs", 0*sizeof(long),  & N_MULTI_WGS, sizeof(N_MULTI_WGS)));
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
						int ss_temp_length = 0;
						if ( ! (Work_Group_Id<Work_Group_Id_Limit))
						{
							long INPUT = -1;
							DPU_ASSERT(dpu_copy_to(ss_dpu, "INPUT", 0,  & INPUT, sizeof(INPUT)));
							continue;
						}
					long long m_dev_map[] = {0};
					long long f_m_dev_stride, m_dev_offset[1] = {-1}, m_dev_end[1] = {-1}, m_dev_size = 0, m_dev_cpu_offset[1] = {-1}, m_dev_length[1] = {0}, m_dev_ele_size;
					long long m_dev_stride[] = {size};
					long long m_dev_flag[] = {2};
					long long a_dev_map[] = {0, 1};
					long long f_a_dev_stride, a_dev_offset[2] = {-1, -1}, a_dev_end[2] = {-1, -1}, a_dev_size = 0, a_dev_cpu_offset[2] = {-1, -1}, a_dev_length[2] = {0, 0}, a_dev_ele_size;
					long long a_dev_stride[] = {size, 1};
					long long a_dev_flag[] = {1, 1};
					
						
						dpu_copy_to(ss_dpu, "_bd", 0,  & _bd, sizeof(_bd));
						DPU_ASSERT(dpu_copy_to(ss_dpu, "m_dev_start", 0,  & dpu_mram_offset, sizeof(dpu_mram_offset)));
						m_dev_ele_size=(((char *)(m+1))-((char *)m));
						if (gxl[0]<((size-1)-t))
						{
						long _gxl[] = {gxl[0], ss_min(gxl[1], ((size-1)-t)-1)};
							m_dev_offset[0]=((size*((_gxl[0]+t)+1))+t);
							m_dev_end[0]=((size*((_gxl[1]+t)+1))+t);
						}
						f_m_dev_stride=size;
						ss_SE = (int(*)[4])generateRectangle(m_dev_offset, m_dev_end, m_dev_stride, m_dev_map, m_dev_flag, f_m_dev_stride, 1);
						
						for (int copy_i = 0; copy_i<1; copy_i ++ )
						{
							m_dev_cpu_offset[copy_i]=m_dev_offset[copy_i];
							if (m_dev_offset[copy_i]<0)
							{
								continue;
							}
							// printf("m_dev_index[%d] : %lld %lld, stride=%lld, SE: %d %d %d %d\n", copy_i, m_dev_offset[copy_i], m_dev_end[copy_i], f_m_dev_stride, ss_SE[copy_i][0], ss_SE[copy_i][1], ss_SE[copy_i][2], ss_SE[copy_i][3]);
							ss_temp_length=(((ss_SE[copy_i][3]-ss_SE[copy_i][1])+1)*m_dev_stride[copy_i]);
							if (m_dev_map[copy_i]==copy_i)
							{
								if (((m_dev_end[copy_i]-m_dev_offset[copy_i])+1)>ss_temp_length)
								{
									int temp_ind = 0;
									float * ss_temp = ((void *)0);
									m_dev_offset[copy_i]=0;
									m_dev_end[copy_i]=(ss_temp_length-1);
									m_dev_flag[copy_i]|=4;
									ss_temp_length=make8ByteAligned(ss_temp_length, sizeof(float));
									if (m_dev_flag[copy_i]&1)
									{
										ss_temp=malloc(ss_temp_length*sizeof(float));
										for (int ss_i = ss_SE[copy_i][1]; ss_i<=ss_SE[copy_i][3]; ss_i+=1)
										{
											for (int ss_j = ss_SE[copy_i][0]; ss_j<=ss_SE[copy_i][2]; ss_j ++ )
											{
												ss_temp[temp_ind ++ ]=m[(size*ss_i)+ss_j];
											}
										}
										dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset, ss_temp, make8ByteAligned((m_dev_end[copy_i]-m_dev_offset[copy_i])+1, m_dev_ele_size)*m_dev_ele_size);
										free(ss_temp);
									}
								}
								else
								{
									int copy_to_size = (make8ByteAligned((m_dev_end[copy_i]-m_dev_offset[copy_i])+1, m_dev_ele_size)*m_dev_ele_size);
									if (m_dev_flag[copy_i]&1)
									{
										dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset, m+m_dev_offset[copy_i], copy_to_size);
									}
									m_dev_stride[copy_i]=1;
								}
								m_dev_length[copy_i]=((m_dev_end[copy_i]-m_dev_offset[copy_i])+1);
								m_dev_offset[copy_i]=m_dev_size;
								m_dev_size=make8ByteAligned(m_dev_length[copy_i], m_dev_ele_size);
							}
							else
							{
								m_dev_stride[copy_i]=m_dev_stride[m_dev_map[copy_i]];
								m_dev_offset[copy_i]=(m_dev_offset[m_dev_map[copy_i]]+((m_dev_stride[copy_i]*ss_SE[copy_i][1])+ss_SE[copy_i][0]));
							}
							// printf("m_dev: %d : %lld %lld L%lld S%lld F%lld M%lld\n", copy_i, m_dev_offset[copy_i], m_dev_cpu_offset[copy_i], m_dev_length[copy_i], m_dev_stride[copy_i], m_dev_flag[copy_i], m_dev_map[copy_i]);
						}
						dpu_copy_to(ss_dpu, "m_dev_offset", 0,  & m_dev_offset, sizeof(m_dev_offset));
						dpu_copy_to(ss_dpu, "m_dev_stride", 0,  & m_dev_stride, sizeof(m_dev_stride));
						dpu_copy_to(ss_dpu, "m_dev_flag", 0,  & m_dev_flag, sizeof(m_dev_flag));
						dpu_copy_to(ss_dpu, "m_dev_length", 0,  & m_dev_length, sizeof(m_dev_length));
						dpu_copy_to(ss_dpu, "m_dev_cpu_offset", 0,  & m_dev_cpu_offset, sizeof(m_dev_cpu_offset));
						dpu_copy_to(ss_dpu, "m_dev_map", 0,  & m_dev_map, sizeof(m_dev_map));
						dpu_mram_offset+=(m_dev_size*m_dev_ele_size);
						
						
						DPU_ASSERT(dpu_copy_to(ss_dpu, "a_dev_start", 0,  & dpu_mram_offset, sizeof(dpu_mram_offset)));
						a_dev_ele_size=(((char *)(a+1))-((char *)a));
						if (gxl[0]<((size-1)-t))
						{
						long _gxl[] = {gxl[0], ss_min(gxl[1], ((size-1)-t)-1)};
							a_dev_offset[0]=((size*((_gxl[0]+t)+1))+t);
							a_dev_end[0]=((size*((_gxl[1]+t)+1))+t);
							a_dev_offset[1]=((size*t)+t);
							a_dev_end[1]=((size*t)+t);
						}
						f_a_dev_stride=size;
						ss_SE = (int(*)[4])generateRectangle(a_dev_offset, a_dev_end, a_dev_stride, a_dev_map, a_dev_flag, f_a_dev_stride, 2);
						
						for (int copy_i = 0; copy_i<2; copy_i ++ )
						{
							a_dev_cpu_offset[copy_i]=a_dev_offset[copy_i];
							if (a_dev_offset[copy_i]<0)
							{
								continue;
							}
							// printf("a_dev_index[%d] : %lld %lld, stride=%lld, SE: %d %d %d %d\n", copy_i, a_dev_offset[copy_i], a_dev_end[copy_i], f_a_dev_stride, ss_SE[copy_i][0], ss_SE[copy_i][1], ss_SE[copy_i][2], ss_SE[copy_i][3]);
							ss_temp_length=(((ss_SE[copy_i][3]-ss_SE[copy_i][1])+1)*a_dev_stride[copy_i]);
							if (a_dev_map[copy_i]==copy_i)
							{
								if (((a_dev_end[copy_i]-a_dev_offset[copy_i])+1)>ss_temp_length)
								{
									int temp_ind = 0;
									float * ss_temp = ((void *)0);
									a_dev_offset[copy_i]=0;
									a_dev_end[copy_i]=(ss_temp_length-1);
									a_dev_flag[copy_i]|=4;
									ss_temp_length=make8ByteAligned(ss_temp_length, sizeof(float));
									if (a_dev_flag[copy_i]&1)
									{
										ss_temp=malloc(ss_temp_length*sizeof(float));
										for (int ss_i = ss_SE[copy_i][1]; ss_i<=ss_SE[copy_i][3]; ss_i+=1)
										{
											for (int ss_j = ss_SE[copy_i][0]; ss_j<=ss_SE[copy_i][2]; ss_j ++ )
											{
												ss_temp[temp_ind ++ ]=a[(size*ss_i)+ss_j];
											}
										}
										dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(a_dev_size*a_dev_ele_size), ss_temp, make8ByteAligned((a_dev_end[copy_i]-a_dev_offset[copy_i])+1, a_dev_ele_size)*a_dev_ele_size);
										free(ss_temp);
									}
								}
								else
								{
									int copy_to_size = (make8ByteAligned((a_dev_end[copy_i]-a_dev_offset[copy_i])+1, a_dev_ele_size)*a_dev_ele_size);
									if (a_dev_flag[copy_i]&1)
									{
										dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(a_dev_size*a_dev_ele_size), a+a_dev_offset[copy_i], copy_to_size);
									}
									a_dev_stride[copy_i]=1;
								}
								a_dev_length[copy_i]=((a_dev_end[copy_i]-a_dev_offset[copy_i])+1);
								a_dev_offset[copy_i]=a_dev_size;
								a_dev_size+=make8ByteAligned(a_dev_length[copy_i], a_dev_ele_size);
							}
							else
							{
								a_dev_stride[copy_i]=a_dev_stride[a_dev_map[copy_i]];
								a_dev_offset[copy_i]=(a_dev_offset[a_dev_map[copy_i]]+((a_dev_stride[copy_i]*ss_SE[copy_i][1])+ss_SE[copy_i][0]));
							}
							// printf("a_dev: %d : %lld %lld L%lld S%lld F%lld M%lld\n", copy_i, a_dev_offset[copy_i], a_dev_cpu_offset[copy_i], a_dev_length[copy_i], a_dev_stride[copy_i], a_dev_flag[copy_i], a_dev_map[copy_i]);
						}
						dpu_copy_to(ss_dpu, "a_dev_offset", 0,  & a_dev_offset, sizeof(a_dev_offset));
						dpu_copy_to(ss_dpu, "a_dev_stride", 0,  & a_dev_stride, sizeof(a_dev_stride));
						dpu_copy_to(ss_dpu, "a_dev_flag", 0,  & a_dev_flag, sizeof(a_dev_flag));
						dpu_copy_to(ss_dpu, "a_dev_length", 0,  & a_dev_length, sizeof(a_dev_length));
						dpu_copy_to(ss_dpu, "a_dev_cpu_offset", 0,  & a_dev_cpu_offset, sizeof(a_dev_cpu_offset));
						dpu_copy_to(ss_dpu, "a_dev_map", 0,  & a_dev_map, sizeof(a_dev_map));
						dpu_mram_offset+=(a_dev_size*a_dev_ele_size);
						
						
						// printf("WG id:%ld bxl: %ld %ld gxl: %ld %ld\n", Work_Group_Id, bxl[0], bxl[1], gxl[0], gxl[1]);
						int n_multi_wgs = bxl[1]-bxl[0]+1;
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
						int ss_temp_length = 0;
						if ( ! (p_Work_Group_Id<Work_Group_Id_Limit))
						{
							break;
						}
						{
							long long f_p_m_dev_stride, p_m_dev_offset[1], p_m_dev_length[1], p_m_dev_stride[1], p_m_dev_cpu_offset[1], p_m_dev_flag[1], p_m_dev_map[1], p_m_dev_ele_size;
							DPU_ASSERT(dpu_copy_from(ss_dpu, "m_dev_start", 0,  & dpu_mram_offset, sizeof(dpu_mram_offset)));
							p_m_dev_ele_size=(((char *)(m+1))-((char *)m));
							dpu_copy_from(ss_dpu, "m_dev_offset", 0,  & p_m_dev_offset, sizeof(p_m_dev_offset));
							dpu_copy_from(ss_dpu, "m_dev_stride", 0,  & p_m_dev_stride, sizeof(p_m_dev_stride));
							dpu_copy_from(ss_dpu, "m_dev_flag", 0,  & p_m_dev_flag, sizeof(p_m_dev_flag));
							dpu_copy_from(ss_dpu, "m_dev_length", 0,  & p_m_dev_length, sizeof(p_m_dev_length));
							dpu_copy_from(ss_dpu, "m_dev_cpu_offset", 0,  & p_m_dev_cpu_offset, sizeof(p_m_dev_cpu_offset));
							dpu_copy_from(ss_dpu, "m_dev_map", 0,  & p_m_dev_map, sizeof(p_m_dev_map));
							f_p_m_dev_stride=size;
							for (int copy_i = 0; copy_i<1; copy_i ++ )
							{
								if (p_m_dev_offset[copy_i]<0)
								{
									continue;
								}
								// printf("p_m_dev: %d : %lld %lld L%lld S%lld F%lld M%lld\n", copy_i, p_m_dev_offset[copy_i], p_m_dev_cpu_offset[copy_i], p_m_dev_length[copy_i], p_m_dev_stride[copy_i], p_m_dev_flag[copy_i], p_m_dev_map[copy_i]);
								ss_temp_length=p_m_dev_length[copy_i];
								if (p_m_dev_map[copy_i]==copy_i)
								{
									if (p_m_dev_flag[copy_i]&4)
									{
										int temp_ind = 0;
										float * ss_temp = ((void *)0);
										ss_temp_length=make8ByteAligned(ss_temp_length, sizeof(float));
										if (p_m_dev_flag[copy_i]&2)
										{
											int ss_i_start = (p_m_dev_cpu_offset[copy_i]/size);
											int ss_i_limit = (ss_i_start+(p_m_dev_length[copy_i]/p_m_dev_stride[copy_i]));
											int ss_j_start = (p_m_dev_cpu_offset[copy_i]%size);
											int ss_j_limit = (p_m_dev_stride[copy_i]+ss_j_start);
											// printf("ss_i:%d-%d, ss_j:%d-%d\n",ss_i_start, ss_i_limit, ss_j_start, ss_j_limit);
											ss_temp=malloc(ss_temp_length*sizeof(float));
											dpu_copy_from(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(p_m_dev_offset[copy_i]*p_m_dev_ele_size), ss_temp, make8ByteAligned(p_m_dev_length[copy_i], p_m_dev_ele_size)*p_m_dev_ele_size);
											for (int ss_i = ss_i_start; ss_i<ss_i_limit; ss_i+=1)
											{
												for (int ss_j = ss_j_start; ss_j<ss_j_limit; ss_j ++ )
												{
													m[(size*ss_i)+ss_j]=ss_temp[temp_ind ++ ];
													if (temp_ind>=p_m_dev_length[copy_i])
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
										int copy_from_size = (make8ByteAligned(p_m_dev_length[copy_i], p_m_dev_ele_size)*p_m_dev_ele_size);
										if (p_m_dev_flag[copy_i]&2)
										{
											float * ss_temp = ((void *)0);
											int temp_ind = 0;
											if (copy_from_size>p_m_dev_length[copy_i])
											{
												ss_temp=malloc((copy_from_size-p_m_dev_length[copy_i])*sizeof(float));
												for (int ss_i = (p_m_dev_cpu_offset[copy_i]+p_m_dev_length[copy_i]); ss_i<copy_from_size; ss_i ++ )
												{
													ss_temp[temp_ind ++ ]=m[ss_i];
												}
											}
											dpu_copy_from(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(p_m_dev_offset[copy_i]*p_m_dev_ele_size), m+p_m_dev_cpu_offset[copy_i], copy_from_size);
											if (copy_from_size>p_m_dev_length[copy_i])
											{
												temp_ind=0;
												for (int ss_i = (p_m_dev_cpu_offset[copy_i]+p_m_dev_length[copy_i]); ss_i<copy_from_size; ss_i ++ )
												{
													m[ss_i]=ss_temp[temp_ind ++ ];
												}
												free(ss_temp);
											}
										}
									}
								}
							}
						}
						int n_multi_wgs = bxl[1]-bxl[0]+1;
						p_bx=((p_bx+n_multi_wgs)%bx_limit);
						p_Work_Group_Id+=n_multi_wgs;
					}
					uint64_t dpu_id = 0;
					int64_t cycles[ss_nr_dpus];
					int64_t barrier_count[ss_nr_dpus][N_TASKLETS];
					int64_t work_done[ss_nr_dpus][N_TASKLETS];
					int64_t T_TD_START[ss_nr_dpus][N_TASKLETS];
					int64_t T_TD_END[ss_nr_dpus][N_TASKLETS];
					int64_t N_WG_ID[ss_nr_dpus][N_TASKLETS];
					int ss_used_ndpus = (Work_Group_Id_Limit+N_MULTI_WGS-1)/N_MULTI_WGS;
					int ss_ndpus = (ss_nr_dpus < ss_used_ndpus) ? ss_nr_dpus : ss_used_ndpus;
					
					int ss_iter_dpu;
					DPU_FOREACH(ss_dpu_set, ss_dpu, ss_iter_dpu) {
						if(ss_iter_dpu>ss_ndpus)
						break;
						// DPU_ASSERT(dpu_log_read(ss_dpu, stdout));
						DPU_ASSERT(dpu_copy_from(ss_dpu, "cycles", 0, &cycles[dpu_id], sizeof(int64_t)));
						DPU_ASSERT(dpu_copy_from(ss_dpu, "barrier_count", 0, &barrier_count[dpu_id], sizeof(barrier_count[0])));
						DPU_ASSERT(dpu_copy_from(ss_dpu, "WORK_DONE", 0, &work_done[dpu_id], sizeof(work_done[0])));
						DPU_ASSERT(dpu_copy_from(ss_dpu, "T_TD_START", 0, &T_TD_START[dpu_id], sizeof(T_TD_START[0])));
						DPU_ASSERT(dpu_copy_from(ss_dpu, "T_TD_END", 0, &T_TD_END[dpu_id], sizeof(T_TD_END[0])));
						DPU_ASSERT(dpu_copy_from(ss_dpu, "N_WG_ID", 0, &N_WG_ID[dpu_id], sizeof(N_WG_ID[0])));
						dpu_id++;
					}
					
					uint64_t total_cycles =0;
					int32_t max_count = 0;
					
					for(int i=0; i<ss_ndpus; i++) {
						// printf("%ld ", cycles[i]);
						total_cycles+=cycles[i];
						if(max_count<cycles[i])
							max_count = cycles[i];
					}
					
					FINAL_TOTAL_CYCLES_0 += max_count;

					printf("Work_Group_Id: %ld t:%ld bx: %ld by: %ld\n", Work_Group_Id, t, bx, by);
					// printf("kernel0 - Total cycles = %" PRId64 "\n", total_cycles);
					printf("kernel0 - Max cycles = %" PRId32 "\n", max_count);
					
					// printf("----------\n");
				int trace_values[5] = {N_TASKLETS, N_MULTI_WGS, PARTITION_DIM_GRID, PARTITION_DIM_WG, max_count};
					note_down(trace_values, 0, t);
					
					for(int i=0; i<ss_ndpus; i++) {
						// printf("Dpu %d:\t",i);
						for(int j=0; j<N_TASKLETS; j++) {
							if(work_done[i][j] == 0) {
								printf("Dpu %d:\t",i);
								printf("[ERROR] WORK NOT DONE\n");
							}
							
						}
						// printf("\n");
					}
					
				}
			}
			
			
			DPU_ASSERT( dpu_free(ss_dpu_set) );
			
		}
		CL_SUCCESS;
		if (timing)
		{
			/*             printf("here1a\n"); */
			kernelTime+=eventTime(kernelEvent, command_queue);
			/*             printf("here1b\n"); */
		}
		CL_SUCCESS;
		/* Fan1<<<dimGrid,dimBlock>>>(m_cuda,a_cuda,Size,t); */
		/* cudaThreadSynchronize(); */
		/* kernel args */
		argchk=CL_SUCCESS;
		argchk|=CL_SUCCESS;
		argchk|=CL_SUCCESS;
		argchk|=CL_SUCCESS;
		argchk|=CL_SUCCESS;
		CL_SUCCESS;
		/* launch kernel */
		{
			
			struct dpu_set_t ss_dpu_set, ss_dpu;
			char* ss_profile = "backend=simulator";DPU_ASSERT(dpu_alloc(-1, ss_profile, &ss_dpu_set));
			uint32_t ss_nr_dpus = 0, ss_nr_ranks=0; 
			dpu_get_nr_dpus(ss_dpu_set, &ss_nr_dpus); 
			
			
			
			
			
			
			long long N_MULTI_WGS = group_size;
			long Work_Group_Id = 0;
			long p_Work_Group_Id = 0;
		long ss_local_work[] = {localWorksizeFan2[0], localWorksizeFan2[1]};
			long bx = 0, p_bx = 0, by = 0, p_by = 0, bz = 0, p_bz = 0;
		long txl[] = {0, ss_local_work[0]};
			long bx_limit = (globalWorksizeFan2[0]/ss_local_work[0]);
		long tyl[] = {0, ss_local_work[1]};
			long by_limit = (globalWorksizeFan2[1]/ss_local_work[1]);
			long Work_Group_Id_Limit = (bx_limit*by_limit);
			long Thread_Id_Limit = (globalWorksizeFan2[0]*globalWorksizeFan2[1]);
			long N_TASKLETS = n_tasklets;
			long PARTITION_DIM_GRID = 0;
			long PARTITION_DIM_WG = 0;
			DPU_ASSERT(dpu_load(ss_dpu_set, "Fan2", NULL));
			
			if(t==size-2 || t==0){
				printf("DPU_ALLOCATED = %d\n", ss_nr_dpus); 
				dpu_get_nr_ranks(ss_dpu_set, &ss_nr_ranks); 
				printf("DPU_RANKS = %d\n", ss_nr_ranks); 
				printf("Work_Group_Id_Limit = %ld\n", Work_Group_Id_Limit);
					printf("Kernel :1 \n");
					printf("Global : ");
					printf("[0]%ld ", globalWorksizeFan2[0]);
					printf("[1]%ld ", globalWorksizeFan2[1]);
					printf("\n");
					printf("Local : ");
					printf("[0]%ld ", ss_local_work[0]);
					printf("[1]%ld ", ss_local_work[1]);
					printf("\n");	
			}
			
			
			DPU_FOREACH(ss_dpu_set, ss_dpu)
			{
				DPU_ASSERT(dpu_copy_to(ss_dpu, "GRID_LIMIT", 0, &bx_limit, sizeof(bx_limit)));
				DPU_ASSERT(dpu_copy_to(ss_dpu, "_gs", 0, globalWorksizeFan2, sizeof(globalWorksizeFan2)));
				DPU_ASSERT(dpu_copy_to(ss_dpu, "_bs", 0, ss_local_work, sizeof(ss_local_work)));
				DPU_ASSERT(dpu_copy_to(ss_dpu, "size", 0,  & size, sizeof(size)));
				DPU_ASSERT(dpu_copy_to(ss_dpu, "t", 0,  & t, sizeof(t)));
				DPU_ASSERT(dpu_copy_to(ss_dpu, "dpu_gs", 0*sizeof(long),  & N_MULTI_WGS, sizeof(N_MULTI_WGS)));
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
					long byl[] = {by, by};
					long gyl[] = {(tyl[1]*byl[0]), ((tyl[1]*byl[1])+(tyl[1]-1))};
					long long _bd[] = {bx, by};
						long long dpu_mram_offset = 0;
						int ss_temp_offset = 0;
						int (*ss_SE)[4];;
						int ss_temp_length = 0;
						if ( ! (Work_Group_Id<Work_Group_Id_Limit))
						{
							long INPUT = -1;
							DPU_ASSERT(dpu_copy_to(ss_dpu, "INPUT", 0,  & INPUT, sizeof(INPUT)));
							continue;
						}
					long long m_dev_map[] = {0, 1};
					long long f_m_dev_stride, m_dev_offset[2] = {-1, -1}, m_dev_end[2] = {-1, -1}, m_dev_size = 0, m_dev_cpu_offset[2] = {-1, -1}, m_dev_length[2] = {0, 0}, m_dev_ele_size;
					long long m_dev_stride[] = {size, size};
					long long m_dev_flag[] = {1, 1};
					long long a_dev_map[] = {0, 1};
					long long f_a_dev_stride, a_dev_offset[2] = {-1, -1}, a_dev_end[2] = {-1, -1}, a_dev_size = 0, a_dev_cpu_offset[2] = {-1, -1}, a_dev_length[2] = {0, 0}, a_dev_ele_size;
					long long a_dev_stride[] = {size, 1};
					long long a_dev_flag[] = {3, 1};
					long long b_dev_map[] = {0, 1};
					long long f_b_dev_stride, b_dev_offset[2] = {-1, -1}, b_dev_end[2] = {-1, -1}, b_dev_size = 0, b_dev_cpu_offset[2] = {-1, -1}, b_dev_length[2] = {0, 0}, b_dev_ele_size;
					long long b_dev_stride[] = {1, 1};
					long long b_dev_flag[] = {3, 1};
					
						dpu_copy_to(ss_dpu, "_bd", 0,  & _bd, sizeof(_bd));
						DPU_ASSERT(dpu_copy_to(ss_dpu, "m_dev_start", 0,  & dpu_mram_offset, sizeof(dpu_mram_offset)));
						m_dev_ele_size=(((char *)(m+1))-((char *)m));
						if ((gxl[0]<((size-1)-t))&&(gyl[0]<(size-t)))
						{
						long _gxl[] = {gxl[0], ss_min(gxl[1], ((size-1)-t)-1)};
						long _gyl[] = {gyl[0], ss_min(gyl[1], (size-t)-1)};
							m_dev_offset[0]=((size*((_gxl[0]+1)+t))+t);
							m_dev_end[0]=((size*((_gxl[1]+1)+t))+t);
							if ((gyl[0]<=0)&&(0<=gyl[1]))
							{
							long _gxl[] = {gxl[0], ss_min(gxl[1], ((size-1)-t)-1)};
							long _gyl[] = {0, ss_min(0, (size-t)-1)};
								m_dev_offset[1]=((size*((_gxl[0]+1)+t))+(_gyl[0]+t));
								m_dev_end[1]=((size*((_gxl[1]+1)+t))+(_gyl[1]+t));
							}
						}
						f_m_dev_stride=size;
						ss_SE = (int(*)[4])generateRectangle(m_dev_offset, m_dev_end, m_dev_stride, m_dev_map, m_dev_flag, f_m_dev_stride, 2);
						
						for (int copy_i = 0; copy_i<2; copy_i ++ )
						{
							m_dev_cpu_offset[copy_i]=m_dev_offset[copy_i];
							if (m_dev_offset[copy_i]<0)
							{
								continue;
							}
							// printf("m_dev_index[%d] : %lld %lld, stride=%lld, SE: %d %d %d %d\n", copy_i, m_dev_offset[copy_i], m_dev_end[copy_i], f_m_dev_stride, ss_SE[copy_i][0], ss_SE[copy_i][1], ss_SE[copy_i][2], ss_SE[copy_i][3]);
							ss_temp_length=(((ss_SE[copy_i][3]-ss_SE[copy_i][1])+1)*m_dev_stride[copy_i]);
							if (m_dev_map[copy_i]==copy_i)
							{
								if (((m_dev_end[copy_i]-m_dev_offset[copy_i])+1)>ss_temp_length)
								{
									int temp_ind = 0;
									float * ss_temp = ((void *)0);
									m_dev_offset[copy_i]=0;
									m_dev_end[copy_i]=(ss_temp_length-1);
									m_dev_flag[copy_i]|=4;
									ss_temp_length=make8ByteAligned(ss_temp_length, sizeof(float));
									if (m_dev_flag[copy_i]&1)
									{
										ss_temp=malloc(ss_temp_length*sizeof(float));
										for (int ss_i = ss_SE[copy_i][1]; ss_i<=ss_SE[copy_i][3]; ss_i+=1)
										{
											for (int ss_j = ss_SE[copy_i][0]; ss_j<=ss_SE[copy_i][2]; ss_j ++ )
											{
												ss_temp[temp_ind ++ ]=m[(size*ss_i)+ss_j];
											}
										}
										dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(m_dev_size*m_dev_ele_size), ss_temp, make8ByteAligned((m_dev_end[copy_i]-m_dev_offset[copy_i])+1, m_dev_ele_size)*m_dev_ele_size);
										free(ss_temp);
									}
								}
								else
								{
									int copy_to_size = (make8ByteAligned((m_dev_end[copy_i]-m_dev_offset[copy_i])+1, m_dev_ele_size)*m_dev_ele_size);
									if (m_dev_flag[copy_i]&1)
									{
										dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(m_dev_size*m_dev_ele_size), m+m_dev_offset[copy_i], copy_to_size);
									}
									m_dev_stride[copy_i]=1;
								}
								m_dev_length[copy_i]=((m_dev_end[copy_i]-m_dev_offset[copy_i])+1);
								m_dev_offset[copy_i]=m_dev_size;
								m_dev_size+=make8ByteAligned(m_dev_length[copy_i], m_dev_ele_size);
							}
							else
							{
								m_dev_stride[copy_i]=m_dev_stride[m_dev_map[copy_i]];
								m_dev_offset[copy_i]=(m_dev_offset[m_dev_map[copy_i]]+((m_dev_stride[copy_i]*ss_SE[copy_i][1])+ss_SE[copy_i][0]));
							}
							// printf("m_dev: %d : %lld %lld L%lld S%lld F%lld M%lld\n", copy_i, m_dev_offset[copy_i], m_dev_cpu_offset[copy_i], m_dev_length[copy_i], m_dev_stride[copy_i], m_dev_flag[copy_i], m_dev_map[copy_i]);
						}
						dpu_copy_to(ss_dpu, "m_dev_offset", 0,  & m_dev_offset, sizeof(m_dev_offset));
						dpu_copy_to(ss_dpu, "m_dev_stride", 0,  & m_dev_stride, sizeof(m_dev_stride));
						dpu_copy_to(ss_dpu, "m_dev_flag", 0,  & m_dev_flag, sizeof(m_dev_flag));
						dpu_copy_to(ss_dpu, "m_dev_length", 0,  & m_dev_length, sizeof(m_dev_length));
						dpu_copy_to(ss_dpu, "m_dev_cpu_offset", 0,  & m_dev_cpu_offset, sizeof(m_dev_cpu_offset));
						dpu_copy_to(ss_dpu, "m_dev_map", 0,  & m_dev_map, sizeof(m_dev_map));
						dpu_mram_offset+=(m_dev_size*m_dev_ele_size);
						
						
						DPU_ASSERT(dpu_copy_to(ss_dpu, "a_dev_start", 0,  & dpu_mram_offset, sizeof(dpu_mram_offset)));
						a_dev_ele_size=(((char *)(a+1))-((char *)a));
						if ((gxl[0]<((size-1)-t))&&(gyl[0]<(size-t)))
						{
						long _gxl[] = {gxl[0], ss_min(gxl[1], ((size-1)-t)-1)};
						long _gyl[] = {gyl[0], ss_min(gyl[1], (size-t)-1)};
							a_dev_offset[0]=((size*((_gxl[0]+1)+t))+(_gyl[0]+t));
							a_dev_end[0]=((size*((_gxl[1]+1)+t))+(_gyl[1]+t));
							a_dev_offset[1]=((size*t)+(_gyl[0]+t));
							a_dev_end[1]=((size*t)+(_gyl[1]+t));
						}
						f_a_dev_stride=size;
						ss_SE = (int(*)[4])generateRectangle(a_dev_offset, a_dev_end, a_dev_stride, a_dev_map, a_dev_flag, f_a_dev_stride, 2);
						
						for (int copy_i = 0; copy_i<2; copy_i ++ )
						{
							a_dev_cpu_offset[copy_i]=a_dev_offset[copy_i];
							if (a_dev_offset[copy_i]<0)
							{
								continue;
							}
							// printf("a_dev_index[%d] : %lld %lld, stride=%lld, SE: %d %d %d %d\n", copy_i, a_dev_offset[copy_i], a_dev_end[copy_i], f_a_dev_stride, ss_SE[copy_i][0], ss_SE[copy_i][1], ss_SE[copy_i][2], ss_SE[copy_i][3]);
							ss_temp_length=(((ss_SE[copy_i][3]-ss_SE[copy_i][1])+1)*a_dev_stride[copy_i]);
							if (a_dev_map[copy_i]==copy_i)
							{
								if (((a_dev_end[copy_i]-a_dev_offset[copy_i])+1)>ss_temp_length)
								{
									int temp_ind = 0;
									float * ss_temp = ((void *)0);
									a_dev_offset[copy_i]=0;
									a_dev_end[copy_i]=(ss_temp_length-1);
									a_dev_flag[copy_i]|=4;
									ss_temp_length=make8ByteAligned(ss_temp_length, sizeof(float));
									if (a_dev_flag[copy_i]&1)
									{
										ss_temp=malloc(ss_temp_length*sizeof(float));
										for (int ss_i = ss_SE[copy_i][1]; ss_i<=ss_SE[copy_i][3]; ss_i+=1)
										{
											for (int ss_j = ss_SE[copy_i][0]; ss_j<=ss_SE[copy_i][2]; ss_j ++ )
											{
												ss_temp[temp_ind ++ ]=a[(size*ss_i)+ss_j];
											}
										}
										dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(a_dev_size*a_dev_ele_size), ss_temp, make8ByteAligned((a_dev_end[copy_i]-a_dev_offset[copy_i])+1, a_dev_ele_size)*a_dev_ele_size);
										free(ss_temp);
									}
								}
								else
								{
									int copy_to_size = (make8ByteAligned((a_dev_end[copy_i]-a_dev_offset[copy_i])+1, a_dev_ele_size)*a_dev_ele_size);
									if (a_dev_flag[copy_i]&1)
									{
										dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(a_dev_size*a_dev_ele_size), a+a_dev_offset[copy_i], copy_to_size);
									}
									a_dev_stride[copy_i]=1;
								}
								a_dev_length[copy_i]=((a_dev_end[copy_i]-a_dev_offset[copy_i])+1);
								a_dev_offset[copy_i]=a_dev_size;
								a_dev_size+=make8ByteAligned(a_dev_length[copy_i], a_dev_ele_size);
							}
							else
							{
								a_dev_stride[copy_i]=a_dev_stride[a_dev_map[copy_i]];
								a_dev_offset[copy_i]=(a_dev_offset[a_dev_map[copy_i]]+((a_dev_stride[copy_i]*ss_SE[copy_i][1])+ss_SE[copy_i][0]));
							}
							// printf("a_dev: %d : %lld %lld L%lld S%lld F%lld M%lld\n", copy_i, a_dev_offset[copy_i], a_dev_cpu_offset[copy_i], a_dev_length[copy_i], a_dev_stride[copy_i], a_dev_flag[copy_i], a_dev_map[copy_i]);
						}
						dpu_copy_to(ss_dpu, "a_dev_offset", 0,  & a_dev_offset, sizeof(a_dev_offset));
						dpu_copy_to(ss_dpu, "a_dev_stride", 0,  & a_dev_stride, sizeof(a_dev_stride));
						dpu_copy_to(ss_dpu, "a_dev_flag", 0,  & a_dev_flag, sizeof(a_dev_flag));
						dpu_copy_to(ss_dpu, "a_dev_length", 0,  & a_dev_length, sizeof(a_dev_length));
						dpu_copy_to(ss_dpu, "a_dev_cpu_offset", 0,  & a_dev_cpu_offset, sizeof(a_dev_cpu_offset));
						dpu_copy_to(ss_dpu, "a_dev_map", 0,  & a_dev_map, sizeof(a_dev_map));
						dpu_mram_offset+=(a_dev_size*a_dev_ele_size);
						
						
						DPU_ASSERT(dpu_copy_to(ss_dpu, "b_dev_start", 0,  & dpu_mram_offset, sizeof(dpu_mram_offset)));
						b_dev_ele_size=(((char *)(b+1))-((char *)b));
						if ((gxl[0]<((size-1)-t))&&(gyl[0]<(size-t)))
						{
							if ((gyl[0]<=0)&&(0<=gyl[1]))
							{
							long _gxl[] = {gxl[0], ss_min(gxl[1], ((size-1)-t)-1)};
							long _gyl[] = {0, 0};
								b_dev_offset[0]=((_gxl[0]+1)+t);
								b_dev_end[0]=((_gxl[1]+1)+t);
								b_dev_offset[1]=t;
								b_dev_end[1]=t;
							}
						}
						f_b_dev_stride=1;
						ss_SE = (int(*)[4])generateRectangle(b_dev_offset, b_dev_end, b_dev_stride, b_dev_map, b_dev_flag, f_b_dev_stride, 2);
						
						for (int copy_i = 0; copy_i<2; copy_i ++ )
						{
							b_dev_cpu_offset[copy_i]=b_dev_offset[copy_i];
							if (b_dev_offset[copy_i]<0)
							{
								continue;
							}
							// printf("b_dev_index[%d] : %lld %lld, stride=%lld, SE: %d %d %d %d\n", copy_i, b_dev_offset[copy_i], b_dev_end[copy_i], f_b_dev_stride, ss_SE[copy_i][0], ss_SE[copy_i][1], ss_SE[copy_i][2], ss_SE[copy_i][3]);
							ss_temp_length=(((ss_SE[copy_i][3]-ss_SE[copy_i][1])+1)*b_dev_stride[copy_i]);
							if (b_dev_map[copy_i]==copy_i)
							{
								if (((b_dev_end[copy_i]-b_dev_offset[copy_i])+1)>ss_temp_length)
								{
									int temp_ind = 0;
									float * ss_temp = ((void *)0);
									b_dev_offset[copy_i]=0;
									b_dev_end[copy_i]=(ss_temp_length-1);
									b_dev_flag[copy_i]|=4;
									ss_temp_length=make8ByteAligned(ss_temp_length, sizeof(float));
									if (b_dev_flag[copy_i]&1)
									{
										ss_temp=malloc(ss_temp_length*sizeof(float));
										for (int ss_i = ss_SE[copy_i][1]; ss_i<=ss_SE[copy_i][3]; ss_i+=1)
										{
											for (int ss_j = ss_SE[copy_i][0]; ss_j<=ss_SE[copy_i][2]; ss_j ++ )
											{
												ss_temp[temp_ind ++ ]=b[(1*ss_i)+ss_j];
											}
										}
										dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(b_dev_size*b_dev_ele_size), ss_temp, make8ByteAligned((b_dev_end[copy_i]-b_dev_offset[copy_i])+1, b_dev_ele_size)*b_dev_ele_size);
										free(ss_temp);
									}
								}
								else
								{
									int copy_to_size = (make8ByteAligned((b_dev_end[copy_i]-b_dev_offset[copy_i])+1, b_dev_ele_size)*b_dev_ele_size);
									if (b_dev_flag[copy_i]&1)
									{
										dpu_copy_to(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(b_dev_size*b_dev_ele_size), b+b_dev_offset[copy_i], copy_to_size);
									}
									b_dev_stride[copy_i]=1;
								}
								b_dev_length[copy_i]=((b_dev_end[copy_i]-b_dev_offset[copy_i])+1);
								b_dev_offset[copy_i]=b_dev_size;
								b_dev_size+=make8ByteAligned(b_dev_length[copy_i], b_dev_ele_size);
							}
							else
							{
								b_dev_stride[copy_i]=b_dev_stride[b_dev_map[copy_i]];
								b_dev_offset[copy_i]=(b_dev_offset[b_dev_map[copy_i]]+((b_dev_stride[copy_i]*ss_SE[copy_i][1])+ss_SE[copy_i][0]));
							}
							// printf("b_dev: %d : %lld %lld L%lld S%lld F%lld M%lld\n", copy_i, b_dev_offset[copy_i], b_dev_cpu_offset[copy_i], b_dev_length[copy_i], b_dev_stride[copy_i], b_dev_flag[copy_i], b_dev_map[copy_i]);
						}
						dpu_copy_to(ss_dpu, "b_dev_offset", 0,  & b_dev_offset, sizeof(b_dev_offset));
						dpu_copy_to(ss_dpu, "b_dev_stride", 0,  & b_dev_stride, sizeof(b_dev_stride));
						dpu_copy_to(ss_dpu, "b_dev_flag", 0,  & b_dev_flag, sizeof(b_dev_flag));
						dpu_copy_to(ss_dpu, "b_dev_length", 0,  & b_dev_length, sizeof(b_dev_length));
						dpu_copy_to(ss_dpu, "b_dev_cpu_offset", 0,  & b_dev_cpu_offset, sizeof(b_dev_cpu_offset));
						dpu_copy_to(ss_dpu, "b_dev_map", 0,  & b_dev_map, sizeof(b_dev_map));
						dpu_mram_offset+=(b_dev_size*b_dev_ele_size);
						
						
						// printf("WG id:%ld bxl: %ld %ld gxl: %ld %ld\n", Work_Group_Id, bxl[0], bxl[1], gxl[0], gxl[1]);
						int n_multi_wgs = bxl[1]-bxl[0]+1;
						bx=((bx+n_multi_wgs)%bx_limit);
						if (bx==0)
						{
							by=((by+1)%by_limit);
						}
						Work_Group_Id+=n_multi_wgs;
						// printf("WG Id = %ld\n", Work_Group_Id);
					}
					dpu_launch(ss_dpu_set, DPU_SYNCHRONOUS);
					DPU_FOREACH(ss_dpu_set, ss_dpu)
					{
					long bxl[] = {p_bx, ss_min(p_bx+(N_MULTI_WGS-1), bx_limit-1)};
					long gxl[] = {(txl[1]*bxl[0]), ((txl[1]*bxl[1])+(txl[1]-1))};
					long byl[] = {p_by, p_by};
					long gyl[] = {(tyl[1]*byl[0]), ((tyl[1]*byl[1])+(tyl[1]-1))};
					long long p_bd[] = {p_bx, p_by};
						long long dpu_mram_offset = 0;
						int ss_temp_offset = 0;
						int (*ss_SE)[4];;
						int ss_temp_length = 0;
						if ( ! (p_Work_Group_Id<Work_Group_Id_Limit))
						{
							break;
						}
						{
							long long f_p_a_dev_stride, p_a_dev_offset[2], p_a_dev_length[2], p_a_dev_stride[2], p_a_dev_cpu_offset[2], p_a_dev_flag[2], p_a_dev_map[2], p_a_dev_ele_size;
							DPU_ASSERT(dpu_copy_from(ss_dpu, "a_dev_start", 0,  & dpu_mram_offset, sizeof(dpu_mram_offset)));
							p_a_dev_ele_size=(((char *)(a+1))-((char *)a));
							dpu_copy_from(ss_dpu, "a_dev_offset", 0,  & p_a_dev_offset, sizeof(p_a_dev_offset));
							dpu_copy_from(ss_dpu, "a_dev_stride", 0,  & p_a_dev_stride, sizeof(p_a_dev_stride));
							dpu_copy_from(ss_dpu, "a_dev_flag", 0,  & p_a_dev_flag, sizeof(p_a_dev_flag));
							dpu_copy_from(ss_dpu, "a_dev_length", 0,  & p_a_dev_length, sizeof(p_a_dev_length));
							dpu_copy_from(ss_dpu, "a_dev_cpu_offset", 0,  & p_a_dev_cpu_offset, sizeof(p_a_dev_cpu_offset));
							dpu_copy_from(ss_dpu, "a_dev_map", 0,  & p_a_dev_map, sizeof(p_a_dev_map));
							f_p_a_dev_stride=size;
							for (int copy_i = 0; copy_i<2; copy_i ++ )
							{
								if (p_a_dev_offset[copy_i]<0)
								{
									continue;
								}
								// printf("p_a_dev: %d : %lld %lld L%lld S%lld F%lld M%lld\n", copy_i, p_a_dev_offset[copy_i], p_a_dev_cpu_offset[copy_i], p_a_dev_length[copy_i], p_a_dev_stride[copy_i], p_a_dev_flag[copy_i], p_a_dev_map[copy_i]);
								ss_temp_length=p_a_dev_length[copy_i];
								if (p_a_dev_map[copy_i]==copy_i)
								{
									if (p_a_dev_flag[copy_i]&4)
									{
										int temp_ind = 0;
										float * ss_temp = ((void *)0);
										ss_temp_length=make8ByteAligned(ss_temp_length, sizeof(float));
										if (p_a_dev_flag[copy_i]&2)
										{
											int ss_i_start = (p_a_dev_cpu_offset[copy_i]/size);
											int ss_i_limit = (ss_i_start+(p_a_dev_length[copy_i]/p_a_dev_stride[copy_i]));
											int ss_j_start = (p_a_dev_cpu_offset[copy_i]%size);
											int ss_j_limit = (p_a_dev_stride[copy_i]+ss_j_start);
											// printf("ss_i:%d-%d, ss_j:%d-%d\n",ss_i_start, ss_i_limit, ss_j_start, ss_j_limit);
											ss_temp=malloc(ss_temp_length*sizeof(float));
											dpu_copy_from(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(p_a_dev_offset[copy_i]*p_a_dev_ele_size), ss_temp, make8ByteAligned(p_a_dev_length[copy_i], p_a_dev_ele_size)*p_a_dev_ele_size);
											for (int ss_i = ss_i_start; ss_i<ss_i_limit; ss_i+=1)
											{
												for (int ss_j = ss_j_start; ss_j<ss_j_limit; ss_j ++ )
												{
													a[(size*ss_i)+ss_j]=ss_temp[temp_ind ++ ];
													if (temp_ind>=p_a_dev_length[copy_i])
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
										int copy_from_size = (make8ByteAligned(p_a_dev_length[copy_i], p_a_dev_ele_size)*p_a_dev_ele_size);
										if (p_a_dev_flag[copy_i]&2)
										{
											float * ss_temp = ((void *)0);
											int temp_ind = 0;
											if (copy_from_size>p_a_dev_length[copy_i])
											{
												ss_temp=malloc((copy_from_size-p_a_dev_length[copy_i])*sizeof(float));
												for (int ss_i = (p_a_dev_cpu_offset[copy_i]+p_a_dev_length[copy_i]); ss_i<copy_from_size; ss_i ++ )
												{
													ss_temp[temp_ind ++ ]=a[ss_i];
												}
											}
											dpu_copy_from(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(p_a_dev_offset[copy_i]*p_a_dev_ele_size), a+p_a_dev_cpu_offset[copy_i], copy_from_size);
											if (copy_from_size>p_a_dev_length[copy_i])
											{
												temp_ind=0;
												for (int ss_i = (p_a_dev_cpu_offset[copy_i]+p_a_dev_length[copy_i]); ss_i<copy_from_size; ss_i ++ )
												{
													a[ss_i]=ss_temp[temp_ind ++ ];
												}
												free(ss_temp);
											}
										}
									}
								}
							}
						}
						{
							long long f_p_b_dev_stride, p_b_dev_offset[2], p_b_dev_length[2], p_b_dev_stride[2], p_b_dev_cpu_offset[2], p_b_dev_flag[2], p_b_dev_map[2], p_b_dev_ele_size;
							DPU_ASSERT(dpu_copy_from(ss_dpu, "b_dev_start", 0,  & dpu_mram_offset, sizeof(dpu_mram_offset)));
							p_b_dev_ele_size=(((char *)(b+1))-((char *)b));
							dpu_copy_from(ss_dpu, "b_dev_offset", 0,  & p_b_dev_offset, sizeof(p_b_dev_offset));
							dpu_copy_from(ss_dpu, "b_dev_stride", 0,  & p_b_dev_stride, sizeof(p_b_dev_stride));
							dpu_copy_from(ss_dpu, "b_dev_flag", 0,  & p_b_dev_flag, sizeof(p_b_dev_flag));
							dpu_copy_from(ss_dpu, "b_dev_length", 0,  & p_b_dev_length, sizeof(p_b_dev_length));
							dpu_copy_from(ss_dpu, "b_dev_cpu_offset", 0,  & p_b_dev_cpu_offset, sizeof(p_b_dev_cpu_offset));
							dpu_copy_from(ss_dpu, "b_dev_map", 0,  & p_b_dev_map, sizeof(p_b_dev_map));
							f_p_b_dev_stride=1;
							for (int copy_i = 0; copy_i<2; copy_i ++ )
							{
								if (p_b_dev_offset[copy_i]<0)
								{
									continue;
								}
								// printf("p_b_dev: %d : %lld %lld L%lld S%lld F%lld M%lld\n", copy_i, p_b_dev_offset[copy_i], p_b_dev_cpu_offset[copy_i], p_b_dev_length[copy_i], p_b_dev_stride[copy_i], p_b_dev_flag[copy_i], p_b_dev_map[copy_i]);
								ss_temp_length=p_b_dev_length[copy_i];
								if (p_b_dev_map[copy_i]==copy_i)
								{
									if (p_b_dev_flag[copy_i]&4)
									{
										int temp_ind = 0;
										float * ss_temp = ((void *)0);
										ss_temp_length=make8ByteAligned(ss_temp_length, sizeof(float));
										if (p_b_dev_flag[copy_i]&2)
										{
											int ss_i_start = (p_b_dev_cpu_offset[copy_i]/1);
											int ss_i_limit = (ss_i_start+(p_b_dev_length[copy_i]/p_b_dev_stride[copy_i]));
											int ss_j_start = (p_b_dev_cpu_offset[copy_i]%1);
											int ss_j_limit = (p_b_dev_stride[copy_i]+ss_j_start);
											// printf("ss_i:%d-%d, ss_j:%d-%d\n",ss_i_start, ss_i_limit, ss_j_start, ss_j_limit);
											ss_temp=malloc(ss_temp_length*sizeof(float));
											dpu_copy_from(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(p_b_dev_offset[copy_i]*p_b_dev_ele_size), ss_temp, make8ByteAligned(p_b_dev_length[copy_i], p_b_dev_ele_size)*p_b_dev_ele_size);
											for (int ss_i = ss_i_start; ss_i<ss_i_limit; ss_i+=1)
											{
												for (int ss_j = ss_j_start; ss_j<ss_j_limit; ss_j ++ )
												{
													b[(1*ss_i)+ss_j]=ss_temp[temp_ind ++ ];
													if (temp_ind>=p_b_dev_length[copy_i])
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
										int copy_from_size = (make8ByteAligned(p_b_dev_length[copy_i], p_b_dev_ele_size)*p_b_dev_ele_size);
										if (p_b_dev_flag[copy_i]&2)
										{
											float * ss_temp = ((void *)0);
											int temp_ind = 0;
											if (copy_from_size>p_b_dev_length[copy_i])
											{
												ss_temp=malloc((copy_from_size-p_b_dev_length[copy_i])*sizeof(float));
												for (int ss_i = (p_b_dev_cpu_offset[copy_i]+p_b_dev_length[copy_i]); ss_i<copy_from_size; ss_i ++ )
												{
													ss_temp[temp_ind ++ ]=b[ss_i];
												}
											}
											dpu_copy_from(ss_dpu, DPU_MRAM_HEAP_POINTER_NAME, dpu_mram_offset+(p_b_dev_offset[copy_i]*p_b_dev_ele_size), b+p_b_dev_cpu_offset[copy_i], copy_from_size);
											if (copy_from_size>p_b_dev_length[copy_i])
											{
												temp_ind=0;
												for (int ss_i = (p_b_dev_cpu_offset[copy_i]+p_b_dev_length[copy_i]); ss_i<copy_from_size; ss_i ++ )
												{
													b[ss_i]=ss_temp[temp_ind ++ ];
												}
												free(ss_temp);
											}
										}
									}
								}
							}
						}
						int n_multi_wgs = bxl[1]-bxl[0]+1;
						p_bx=((p_bx+n_multi_wgs)%bx_limit);
						if (p_bx==0)
						{
							p_by=((p_by+1)%by_limit);
						}
						p_Work_Group_Id+=n_multi_wgs;
					}
					uint64_t dpu_id = 0;
					int64_t cycles[ss_nr_dpus];
					int64_t barrier_count[ss_nr_dpus][N_TASKLETS];
					int64_t work_done[ss_nr_dpus][N_TASKLETS];
					int64_t T_TD_START[ss_nr_dpus][N_TASKLETS];
					int64_t T_TD_END[ss_nr_dpus][N_TASKLETS];
					int64_t N_WG_ID[ss_nr_dpus][N_TASKLETS];
					int ss_used_ndpus = (Work_Group_Id_Limit+N_MULTI_WGS-1)/N_MULTI_WGS;
					int ss_ndpus = (ss_nr_dpus < ss_used_ndpus) ? ss_nr_dpus : ss_used_ndpus;
					
					int ss_iter_dpu;
					DPU_FOREACH(ss_dpu_set, ss_dpu, ss_iter_dpu) {
						if(ss_iter_dpu>ss_ndpus)
						break;
						// DPU_ASSERT(dpu_log_read(ss_dpu, stdout));
						DPU_ASSERT(dpu_copy_from(ss_dpu, "cycles", 0, &cycles[dpu_id], sizeof(int64_t)));
						DPU_ASSERT(dpu_copy_from(ss_dpu, "barrier_count", 0, &barrier_count[dpu_id], sizeof(barrier_count[0])));
						DPU_ASSERT(dpu_copy_from(ss_dpu, "WORK_DONE", 0, &work_done[dpu_id], sizeof(work_done[0])));
						DPU_ASSERT(dpu_copy_from(ss_dpu, "T_TD_START", 0, &T_TD_START[dpu_id], sizeof(T_TD_START[0])));
						DPU_ASSERT(dpu_copy_from(ss_dpu, "T_TD_END", 0, &T_TD_END[dpu_id], sizeof(T_TD_END[0])));
						DPU_ASSERT(dpu_copy_from(ss_dpu, "N_WG_ID", 0, &N_WG_ID[dpu_id], sizeof(N_WG_ID[0])));
						dpu_id++;
					}
					
					uint64_t total_cycles =0;
					int32_t max_count = 0;
					
					for(int i=0; i<ss_ndpus; i++) {
						// printf("%ld ", cycles[i]);
						total_cycles+=cycles[i];
						if(max_count<cycles[i])
							max_count = cycles[i];
					}
					
					printf("Work_Group_Id: %ld\n", Work_Group_Id);
					// printf("kernel1- Total cycles = %" PRId64 "\n", total_cycles);
					printf("kernel1 - Max cycles = %" PRId32 "\n", max_count);
					
					printf("----------\n");
					FINAL_TOTAL_CYCLES_1 += max_count;

				int trace_values[5] = {N_TASKLETS, N_MULTI_WGS, PARTITION_DIM_GRID, PARTITION_DIM_WG, max_count};
					note_down(trace_values, 1, t);
					
					for(int i=0; i<ss_ndpus; i++) {
						
						for(int j=0; j<N_TASKLETS; j++) {
							if(work_done[i][j] == 0) {
								printf("Dpu %d:\t",i);
								printf("[ERROR] WORK NOT DONE\n");
							}
									
						}
						// printf("\n");
					}
					
				}
			}
			
			
			DPU_ASSERT( dpu_free(ss_dpu_set) );
			
		}
		CL_SUCCESS;
		if (timing)
		{
			/*             printf("here2a\n"); */
			kernelTime+=eventTime(kernelEvent, command_queue);
			/*             printf("here2b\n"); */
		}
		CL_SUCCESS;
		/* Fan2<<<dimGridXY,dimBlockXY>>>(m_cuda,a_cuda,b_cuda,Size,Size-t,t); */
		/* cudaThreadSynchronize(); */
	}

	printf("FINAL_TOTAL_CYCLES=%lld %lld\n", FINAL_TOTAL_CYCLES_0, FINAL_TOTAL_CYCLES_1);
	fprintf(filePointer1, "%d %d %lld\n", group_size, n_tasklets, FINAL_TOTAL_CYCLES_0);
	fprintf(filePointer2, "%d %d %lld\n", group_size, n_tasklets, FINAL_TOTAL_CYCLES_1);
	fprintf(filePointer3, "%d %d %lld\n", group_size, n_tasklets, FINAL_TOTAL_CYCLES_0+FINAL_TOTAL_CYCLES_1);
	fclose(filePointer1);
	fclose(filePointer2);
	fclose(filePointer3);
	/* 5. transfer data off of device */
	/* change to 0 for nonblocking write */
	/* offset */
	error=CL_SUCCESS;
	CL_SUCCESS;
	if (timing)
	{
		readTime+=eventTime(readEvent, command_queue);
	}
	CL_SUCCESS;
	/* change to 0 for nonblocking write */
	/* offset */
	error=CL_SUCCESS;
	CL_SUCCESS;
	if (timing)
	{
		readTime+=eventTime(readEvent, command_queue);
	}
	CL_SUCCESS;
	/* change to 0 for nonblocking write */
	/* offset */
	error=CL_SUCCESS;
	CL_SUCCESS;
	if (timing)
	{
		readTime+=eventTime(readEvent, command_queue);
	}
	CL_SUCCESS;
	readMB=((float)(((sizeof (float)*size)*((size+size)+1))/1000000.0));
	return ;
}

float eventTime(cl_event event, cl_command_queue command_queue)
{
	cl_int error = 0;
	cl_ulong eventStart, eventEnd;
	float _ret_val_0;
	CL_SUCCESS;
	error=CL_SUCCESS;
	CL_SUCCESS;
	error=CL_SUCCESS;
	CL_SUCCESS;
	_ret_val_0=((float)((eventEnd-eventStart)/1.0E9));
	return _ret_val_0;
}

/* Ke Wang add a function to generate input internally */
int parseCommandline(int argc, char * argv[], char * filename, int * q, int * t, int * p, int * d, int * size)
{
	int i;
	char flag;
	int _ret_val_0;
	if (argc<2)
	{
		_ret_val_0=1;
		return _ret_val_0;
	}
	/* error */
	/* strncpy(filename,argv[1],100); */
	#pragma loop name parseCommandline#0 
	for (i=1; i<argc; i ++ )
	{
		if (argv[i][0]=='-')
		{
			/* flag */
			flag=argv[i][1];
			switch (flag)
			{
				/* platform */
				case 's':
				i ++ ;
				( * size)=atoi(argv[i]);
				printf("Create matrix internally in parse, size = %d \n",  * size);
				break;
				/* platform */
				case 'f':
				i ++ ;
				strncpy(filename, argv[i], 100);
				printf("Read file from %s \n", filename);
				break;
				/* help */
				case 'h':
				_ret_val_0=1;
				return _ret_val_0;
				break;
				/* quiet */
				case 'q':
				( * q)=1;
				break;
				/* timing */
				case 't':
				( * t)=1;
				break;
				/* platform */
				case 'p':
				i ++ ;
				( * p)=atoi(argv[i]);
				break;
				/* device */
				case 'd':
				i ++ ;
				( * d)=atoi(argv[i]);
				break;
				case 'n':
				i++;
		                n_tasklets = atoi(argv[i]);
		                printf("n_tasklets=%d", n_tasklets);
		                break;
		            	case 'g':
		            	i++;
		                group_size = atoi(argv[i]);
		                printf("group_size=%d", group_size);
		                break;
			}
		}
	}
	/* both p and d must be specified if either are specified */
	if (((( * d)>=0)&&(( * p)<0))||((( * p)>=0)&&(( * d)<0)))
	{
		_ret_val_0=1;
		return _ret_val_0;
	}
	_ret_val_0=0;
	return _ret_val_0;
}

void printUsage()
{
	printf("Gaussian Elimination Usage\n");
	printf("\n");
	printf("gaussianElimination [filename] [-hqt] [-p [int] -d [int]]\n");
	printf("\n");
	printf("example:\n");
	printf("$ ./gaussianElimination matrix4.txt\n");
	printf("\n");
	printf("filename     the filename that holds the matrix data\n");
	printf("\n");
	printf("-h           Display the help file\n");
	printf("-q           Quiet mode. Suppress all text output.\n");
	printf("-t           Print timing information.\n");
	printf("\n");
	printf("-p [int]     Choose the platform (must choose both platform and device)\n");
	printf("-d [int]     Choose the device (must choose both platform and device)\n");
	printf("\n");
	printf("\n");
	printf("Notes: 1. The filename is required as the first parameter.\n");
	printf("       2. If you declare either the device or the platform,\n");
	printf("          you must declare both.\n\n");
	return ;
}

/*
------------------------------------------------------

 InitPerRun() -- Initialize the contents of the

 ** multipier matrix **m

 **------------------------------------------------------


*/
void InitPerRun(int size, float * m)
{
	int i;
	#pragma loop name InitPerRun#0 
	for (i=0; i<(size*size); i ++ )
	{
		( * (m+i))=0.0;
	}
	return ;
}

void BackSub(float * a, float * b, float * finalVec, int size)
{
	/* solve "bottom up" */
	int i, j;
	#pragma loop name BackSub#0 
	for (i=0; i<size; i ++ )
	{
		finalVec[(size-i)-1]=b[(size-i)-1];
		#pragma loop name BackSub#0#0 
		for (j=0; j<i; j ++ )
		{
			finalVec[(size-i)-1]-=(( * ((a+(size*((size-i)-1)))+((size-j)-1)))*finalVec[(size-j)-1]);
		}
		finalVec[(size-i)-1]=(finalVec[(size-i)-1]/( * ((a+(size*((size-i)-1)))+((size-i)-1))));
	}
	return ;
}

void InitMat(FILE * fp, int size, float * ary, int nrow, int ncol)
{
	int i, j;
	#pragma loop name InitMat#0 
	for (i=0; i<nrow; i ++ )
	{
		#pragma loop name InitMat#0#0 
		for (j=0; j<ncol; j ++ )
		{
			fscanf(fp, "%f", (ary+(size*i))+j);
		}
	}
	return ;
}

/*
------------------------------------------------------

 InitAry() -- Initialize the array (vector) by reading

 ** data from the data file

 **------------------------------------------------------


*/
void InitAry(FILE * fp, float * ary, int ary_size)
{
	int i;
	#pragma loop name InitAry#0 
	for (i=0; i<ary_size; i ++ )
	{
		fscanf(fp, "%f",  & ary[i]);
	}
	return ;
}

/*
------------------------------------------------------

 PrintMat() -- Print the contents of the matrix

 **------------------------------------------------------


*/
void PrintMat(float * ary, int size, int nrow, int ncol)
{
	int i, j;
	#pragma loop name PrintMat#0 
	for (i=0; i<nrow; i ++ )
	{
		#pragma loop name PrintMat#0#0 
		for (j=0; j<ncol; j ++ )
		{
			printf("%8.2e ",  * ((ary+(size*i))+j));
		}
		printf("\n");
	}
	printf("\n");
	return ;
}

/*
------------------------------------------------------

 PrintAry() -- Print the contents of the array (vector)

 **------------------------------------------------------


*/
void PrintAry(float * ary, int ary_size)
{
	int i;
	#pragma loop name PrintAry#0 
	for (i=0; i<ary_size; i ++ )
	{
		printf("%.2e ", ary[i]);
	}
	printf("\n\n");
	return ;
}
