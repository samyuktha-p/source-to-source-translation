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
 * limited to the U.S. Export Administration Regulations (“EAR”), (15 C.F.R.  *
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
 * of Industry and Security’s website at http://www.bis.doc.gov/.             *
 \****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//#include <math.h>

//#include <CL/cl.h>

#include "clutils.h"
#include "utils.h"


// The following variables have file scope to simplify
// the utility functions

//! All discoverable OpenCL platforms
static cl_platform_id* platforms = NULL;
static cl_uint numPlatforms;

//! All discoverable OpenCL devices (one pointer per platform)
static cl_device_id** devices = NULL;
static cl_uint* numDevices;

//! The chosen OpenCL platform
static cl_platform_id platform = NULL;

//! The chosen OpenCL device
static cl_device_id device = NULL;

//! OpenCL context
static cl_context context = NULL;        

//! OpenCL command queue
static cl_command_queue commandQueue = NULL;  
static cl_command_queue commandQueueProf = NULL;
static cl_command_queue commandQueueNoProf = NULL;

//! Global status of events
static short eventsEnabled = false;

//-------------------------------------------------------
//          Utility functions
//-------------------------------------------------------

//! Take a string and an int, and return a string
char* catStringWithInt(const char* string, int integer) {
    
    if(integer > 99999) {
        printf("Can't handle event identifiers with 6 digits\n");
        exit(-1);
    }

    // 5 characters for the identifier, 1 for the null terminator
    int strLen = strlen(string)+5+1;
    char* eventStr = (char*)alloc(sizeof(char)*strLen);

    char tmp[6];

    strcpy(eventStr, string);
    strncat(eventStr, itoa_portable(integer, tmp, 10), 5);

    return eventStr;
}

/**
 ** C++ version 0.4 char* style "itoa":
 ** Written by Lukás Chmela
 ** Released under GPLv3.
 **/
//portable itoa function
char* itoa_portable(int value, char* result, int base) {
    // check that the base if valid
    if (base < 2 || base > 36) { *result = '\0'; return result; }

    char* ptr = result, *ptr1 = result, tmp_char;
    int tmp_value;

    do {
        tmp_value = value;
        value /= base;
        *ptr++ = "zyxwvutsrqponmlkjihgfedcba9876543210123456789abcdefghijklmnopqrstuvwxyz" [35 + (tmp_value - value * base)];
    } while ( value );

    //Apply negative sign
    if (tmp_value < 0) *ptr++ = '-';
    *ptr-- = '\0';

    while(ptr1 < ptr) {
        tmp_char = *ptr;
        *ptr--= *ptr1;
        *ptr1++ = tmp_char;
    }

    return result;
}

//-------------------------------------------------------
//          Error handling
//-------------------------------------------------------

//! OpenCl error code list
/*!
    An array of character strings used to give the error corresponding to the error code \n

    The error code is the index within this array
*/
char *cl_errs[MAX_ERR_VAL] = {
    (char *)"CL_SUCCESS",                         // 0                            
    (char *)"CL_DEVICE_NOT_FOUND",                //-1                         
    (char *)"CL_DEVICE_NOT_AVAILABLE",            //-2                    
    (char *)"CL_COMPILER_NOT_AVAILABLE",          //-3                 
    (char *)"CL_MEM_OBJECT_ALLOCATION_FAILURE",   //-4            
    (char *)"CL_OUT_OF_RESOURCES",                //-5                         
    (char *)"CL_OUT_OF_HOST_MEMORY",              //-6                      
    (char *)"CL_PROFILING_INFO_NOT_AVAILABLE",    //-7            
    (char *)"CL_MEM_COPY_OVERLAP",                //-8                        
    (char *)"CL_IMAGE_FORMAT_MISMATCH",           //-9                   
    (char *)"CL_IMAGE_FORMAT_NOT_SUPPORTED",      //-10
    (char *)"CL_BUILD_PROGRAM_FAILURE",           //-11           
    (char *)"CL_MAP_FAILURE",                     //-12
    (char *)"",                                   //-13
    (char *)"",                                   //-14
    (char *)"",                                   //-15
    (char *)"",                                   //-16
    (char *)"",                                   //-17
    (char *)"",                                   //-18
    (char *)"",                                   //-19
    (char *)"",                                   //-20
    (char *)"",                                   //-21
    (char *)"",                                   //-22
    (char *)"",                                   //-23
    (char *)"",                                   //-24
    (char *)"",                                   //-25
    (char *)"",                                   //-26
    (char *)"",                                   //-27
    (char *)"",                                   //-28
    (char *)"",                                   //-29
    (char *)"CL_INVALID_VALUE",                   //-30
    (char *)"CL_INVALID_DEVICE_TYPE",             //-31
    (char *)"CL_INVALID_PLATFORM",                //-32
    (char *)"CL_INVALID_DEVICE",                  //-33
    (char *)"CL_INVALID_CONTEXT",                 //-34
    (char *)"CL_INVALID_QUEUE_PROPERTIES",        //-35
    (char *)"CL_INVALID_COMMAND_QUEUE",           //-36
    (char *)"CL_INVALID_HOST_PTR",                //-37
    (char *)"CL_INVALID_MEM_OBJECT",              //-38
    (char *)"CL_INVALID_IMAGE_FORMAT_DESCRIPTOR", //-39
    (char *)"CL_INVALID_IMAGE_SIZE",              //-40
    (char *)"CL_INVALID_SAMPLER",                 //-41
    (char *)"CL_INVALID_BINARY",                  //-42
    (char *)"CL_INVALID_BUILD_OPTIONS",           //-43
    (char *)"CL_INVALID_PROGRAM",                 //-44
    (char *)"CL_INVALID_PROGRAM_EXECUTABLE",      //-45
    (char *)"CL_INVALID_KERNEL_NAME",             //-46
    (char *)"CL_INVALID_KERNEL_DEFINITION",       //-47
    (char *)"CL_INVALID_KERNEL",                  //-48
    (char *)"CL_INVALID_ARG_INDEX",               //-49
    (char *)"CL_INVALID_ARG_VALUE",               //-50
    (char *)"CL_INVALID_ARG_SIZE",                //-51
    (char *)"CL_INVALID_KERNEL_ARGS",             //-52
    (char *)"CL_INVALID_WORK_DIMENSION ",         //-53
    (char *)"CL_INVALID_WORK_GROUP_SIZE",         //-54
    (char *)"CL_INVALID_WORK_ITEM_SIZE",          //-55
    (char *)"CL_INVALID_GLOBAL_OFFSET",           //-56
    (char *)"CL_INVALID_EVENT_WAIT_LIST",         //-57
    (char *)"CL_INVALID_EVENT",                   //-58
    (char *)"CL_INVALID_OPERATION",               //-59
    (char *)"CL_INVALID_GL_OBJECT",               //-60
    (char *)"CL_INVALID_BUFFER_SIZE",             //-61
    (char *)"CL_INVALID_MIP_LEVEL",               //-62
    (char *)"CL_INVALID_GLOBAL_WORK_SIZE"};       //-63

//! OpenCl Error checker
/*!
Checks for error code as per cl_int returned by OpenCl
\param status Error value as cl_int
\param msg User provided error message 
\return True if Error Seen, False if no error
*/
