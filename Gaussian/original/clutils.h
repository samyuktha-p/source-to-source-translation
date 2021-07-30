/****************************************************************************\ 
 * Copyright (c) 2011, Advanced Micro Devices, Inc.                           *
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
 * limited to the U.S. Export Administration Regulations (EAR), (15 C.F.R.  *
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
 * of Industry and Securitys website at http://www.bis.doc.gov/.             *
 \****************************************************************************/

#ifndef __CL_UTILS_H__
#define __CL_UTILS_H__

//#include <CL/cl.h>

#include "open_cl.h"
// The cl_time type is OS specific
/*#ifdef _WIN32
#include <tchar.h>
#include <Windows.h>
typedef __int64 cl_time; 
#else*/
#include <sys/time.h>
typedef double cl_time;
//#endif


//-------------------------------------------------------
// Utility functions
//-------------------------------------------------------

char* catStringWithInt(const char* str, int integer);

char* itoa_portable(int value, char* result, int base);

//-------------------------------------------------------
// Data types
//-------------------------------------------------------
typedef struct{
    int x;
    int y;
} int2;

typedef struct{
    float x;
    float y;
}float2;

typedef struct{
    float x;
    float y;
    float z;
    float w;
}float4;

//-------------------------------------------------------
// Defines
//-------------------------------------------------------

#define MAX_ERR_VAL 64

#define NUM_PROGRAMS 7

#define NUM_KERNELS 13
#define KERNEL_INIT_DET 0 
#define KERNEL_BUILD_DET 1 
#define KERNEL_SURF_DESC 2
#define KERNEL_NORM_DESC 3
#define KERNEL_NON_MAX_SUP 4
#define KERNEL_GET_ORIENT1 5
#define KERNEL_GET_ORIENT2 6
#define KERNEL_NN 7
#define KERNEL_SCAN 8
#define KERNEL_SCAN4 9
#define KERNEL_TRANSPOSE 10
#define KERNEL_SCANIMAGE 11
#define KERNEL_TRANSPOSEIMAGE 12

#endif
