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
/*
 * Copyright (c) 2011, Advanced Micro Devices, Inc.                           *
 *
 * All rights reserved.                                                       *
 *                                                                            *
 * Redistribution and use in source and binary forms, with or without         *
 * modification, are permitted provided that the following conditions         *
 * are met:                                                                   *
 *                                                                            *
 * Redistributions of source code must retain the above copyright notice,     *
 * this list of conditions and the following disclaimer.                      *
 *                                                                            *
 * Redistributions in binary form must reproduce the above copyright notice,  *
 * this list of conditions and the following disclaimer in the documentation  *
 * and/or other materials provided with the distribution.                     *
 *                                                                            *
 * Neither the name of the copyright holder nor the names of its contributors *
 * may be used to endorse or promote products derived from this software      *
 * without specific prior written permission.                                 *
 *                                                                            *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS        *
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED  *
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR *
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR          *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,      *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,        *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR         *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF     *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING       *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS         *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.               *
 *                                                                            *
 * If you use the software (in whole or in part), you shall adhere to all     *
 * applicable U.S., European, and other export laws, including but not        *
 * limited to the U.S. Export Administration Regulations (Ã¢ÂÂEARÃ¢ÂÂ), (15 C.F.R.  *
 * Sections 730 through 774), and E.U. Council Regulation (EC) No 1334/2000   *
 * of 22 June 2000.  Further, pursuant to Section 740.6 of the EAR, you       *
 * hereby certify that, except pursuant to a license granted by the United    *
 * States Department of Commerce Bureau of Industry and Security or as        *
 * otherwise permitted pursuant to a License Exception under the U.S. Export  *
 * Administration Regulations ("EAR"), you will not (1) export, re-export or  *
 * release to a national of a country in Country Groups D:1, E:1 or E:2 any   *
 * restricted technology, software, or source code you receive hereunder,     *
 * or (2) export to Country Groups D:1, E:1 or E:2 the direct product of such *
 * technology or software, if such foreign produced direct product is subject *
 * to national security controls as identified on the Commerce Control List   *
 *(currently found in Supplement 1 to Part 774 of EAR).  For the most current *
 * Country Group listings, or for additional information about the EAR or     *
 * your obligations under those regulations, please refer to the U.S. Bureau  *
 * of Industry and SecurityÃ¢ÂÂs website at http:www.bis.doc.gov/.             *
 \
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/* #include <math.h> */
/* #include <CLcl.h> */
#include "clutils.h"
#include "utils.h"
/* The following variables have file scope to simplify */
/* the utility functions */
/* ! All discoverable OpenCL platforms */
static cl_platform_id * platforms = (void * )0;
static cl_uint numPlatforms;
/* ! All discoverable OpenCL devices (one pointer per platform) */
static cl_device_id * * devices = (void * )0;
static cl_uint * numDevices;
/* ! The chosen OpenCL platform */
static cl_platform_id platform = (void * )0;
/* ! The chosen OpenCL device */
static cl_device_id device = (void * )0;
/* ! OpenCL context */
static cl_context context = (void * )0;
/* ! OpenCL command queue */
static cl_command_queue commandQueue = (void * )0;
static cl_command_queue commandQueueProf = (void * )0;
static cl_command_queue commandQueueNoProf = (void * )0;
/* ! Global status of events */
static short eventsEnabled = 0;
/* ------------------------------------------------------- */
/*          Utility functions */
/* ------------------------------------------------------- */
/* ! Take a string and an int, and return a string */
char *catStringWithInt(const char * string, int integer)
{
	int strLen = (strlen(string)+5)+1;
	char * eventStr = (char * )alloc(sizeof (char)*strLen);
	char tmp[6];
	if (integer>99999)
	{
		printf("Can't handle event identifiers with 6 digits\n");
		exit( - 1);
	}
	/* 5 characters for the identifier, 1 for the null terminator */
	strcpy(eventStr, string);
	strncat(eventStr, itoa_portable(integer, tmp, 10), 5);
	return eventStr;
}

/*

 ** C++ version 0.4 char* style "itoa":
 ** Written by LukÃÂ¡s Chmela
 ** Released under GPLv3.

*/
/* portable itoa function */
char *itoa_portable(int value, char * result, int base)
{
	/* check that the base if valid */
	char * ptr = result, * ptr1 = result, tmp_char;
	int tmp_value;
	if ((base<2)||(base>36))
	{
		( * result)='\0';
		return result;
	}
	do
	{
		tmp_value=value;
		value/=base;
		( * (ptr ++ ))="zyxwvutsrqponmlkjihgfedcba9876543210123456789abcdefghijklmnopqrstuvwxyz"[35+(tmp_value-(value*base))];
	}while(value);
	
	/* Apply negative sign */
	if (tmp_value<0)
	{
		( * (ptr ++ ))='-';
	}
	( * (ptr -- ))='\0';
	while (ptr1<ptr)
	{
		tmp_char=( * ptr);
		( * (ptr -- ))=( * ptr1);
		( * (ptr1 ++ ))=tmp_char;
	}
	return result;
}

/* ------------------------------------------------------- */
/*          Error handling */
/* ------------------------------------------------------- */
/* ! OpenCl error code list */
/*
!
    An array of character strings used to give the error corresponding to the error code \n

    The error code is the index within this array

*/
/* 0                             */
/* -1                          */
/* -2                     */
/* -3                  */
/* -4             */
/* -5                          */
/* -6                       */
/* -7             */
/* -8                         */
/* -9                    */
/* -10 */
/* -11            */
/* -12 */
/* -13 */
/* -14 */
/* -15 */
/* -16 */
/* -17 */
/* -18 */
/* -19 */
/* -20 */
/* -21 */
/* -22 */
/* -23 */
/* -24 */
/* -25 */
/* -26 */
/* -27 */
/* -28 */
/* -29 */
/* -30 */
/* -31 */
/* -32 */
/* -33 */
/* -34 */
/* -35 */
/* -36 */
/* -37 */
/* -38 */
/* -39 */
/* -40 */
/* -41 */
/* -42 */
/* -43 */
/* -44 */
/* -45 */
/* -46 */
/* -47 */
/* -48 */
/* -49 */
/* -50 */
/* -51 */
/* -52 */
/* -53 */
/* -54 */
/* -55 */
/* -56 */
/* -57 */
/* -58 */
/* -59 */
/* -60 */
/* -61 */
/* -62 */
char * cl_errs[64] = {(char * )"CL_SUCCESS", (char * )"CL_DEVICE_NOT_FOUND", (char * )"CL_DEVICE_NOT_AVAILABLE", (char * )"CL_COMPILER_NOT_AVAILABLE", (char * )"CL_MEM_OBJECT_ALLOCATION_FAILURE", (char * )"CL_OUT_OF_RESOURCES", (char * )"CL_OUT_OF_HOST_MEMORY", (char * )"CL_PROFILING_INFO_NOT_AVAILABLE", (char * )"CL_MEM_COPY_OVERLAP", (char * )"CL_IMAGE_FORMAT_MISMATCH", (char * )"CL_IMAGE_FORMAT_NOT_SUPPORTED", (char * )"CL_BUILD_PROGRAM_FAILURE", (char * )"CL_MAP_FAILURE", (char * )"", (char * )"", (char * )"", (char * )"", (char * )"", (char * )"", (char * )"", (char * )"", (char * )"", (char * )"", (char * )"", (char * )"", (char * )"", (char * )"", (char * )"", (char * )"", (char * )"", (char * )"CL_INVALID_VALUE", (char * )"CL_INVALID_DEVICE_TYPE", (char * )"CL_INVALID_PLATFORM", (char * )"CL_INVALID_DEVICE", (char * )"CL_INVALID_CONTEXT", (char * )"CL_INVALID_QUEUE_PROPERTIES", (char * )"CL_INVALID_COMMAND_QUEUE", (char * )"CL_INVALID_HOST_PTR", (char * )"CL_INVALID_MEM_OBJECT", (char * )"CL_INVALID_IMAGE_FORMAT_DESCRIPTOR", (char * )"CL_INVALID_IMAGE_SIZE", (char * )"CL_INVALID_SAMPLER", (char * )"CL_INVALID_BINARY", (char * )"CL_INVALID_BUILD_OPTIONS", (char * )"CL_INVALID_PROGRAM", (char * )"CL_INVALID_PROGRAM_EXECUTABLE", (char * )"CL_INVALID_KERNEL_NAME", (char * )"CL_INVALID_KERNEL_DEFINITION", (char * )"CL_INVALID_KERNEL", (char * )"CL_INVALID_ARG_INDEX", (char * )"CL_INVALID_ARG_VALUE", (char * )"CL_INVALID_ARG_SIZE", (char * )"CL_INVALID_KERNEL_ARGS", (char * )"CL_INVALID_WORK_DIMENSION ", (char * )"CL_INVALID_WORK_GROUP_SIZE", (char * )"CL_INVALID_WORK_ITEM_SIZE", (char * )"CL_INVALID_GLOBAL_OFFSET", (char * )"CL_INVALID_EVENT_WAIT_LIST", (char * )"CL_INVALID_EVENT", (char * )"CL_INVALID_OPERATION", (char * )"CL_INVALID_GL_OBJECT", (char * )"CL_INVALID_BUFFER_SIZE", (char * )"CL_INVALID_MIP_LEVEL", (char * )"CL_INVALID_GLOBAL_WORK_SIZE"};
