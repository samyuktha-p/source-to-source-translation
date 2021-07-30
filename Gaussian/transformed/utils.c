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
 * limited to the U.S. Export Administration Regulations (ÃÂEARÃÂ), (15 C.F.R.  *
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
 * of Industry and SecurityÃÂs website at http:www.bis.doc.gov/.             *
 \
*/
#include <stdio.h>
#include <sys/stat.h>
#include <string.h>
#include <stdlib.h>
#include "utils.h"
static short usingImages = 0;
/* ! A wrapper for malloc that checks the return value */
void *alloc(size_t size)
{
	void * ptr = (void * )0;
	ptr=malloc(size);
	if (ptr==((void * )0))
	{
		perror("malloc");
		exit( - 1);
	}
	return ptr;
}

/* This function checks to make sure a file exists before we open it */
void checkFile(char * filename)
{
	struct stat fileStatus;
	if (stat(filename,  & fileStatus)!=0)
	{
		printf("Error opening file: %s\n", filename);
		exit( - 1);
	}
	else
	{
		if ( ! (32768&fileStatus.st_mode))
		{
			printf("File %s is not a regular file\n", filename);
			exit( - 1);
		}
	}
	return ;
}

/* This function checks to make sure a directory exists  */
void checkDir(char * dirpath)
{
	struct stat fileStatus;
	if (stat(dirpath,  & fileStatus)!=0)
	{
		printf("Directory does not exist: %s\n", dirpath);
		exit( - 1);
	}
	else
	{
		if ( ! (16384&fileStatus.st_mode))
		{
			printf("Directory was not provided: %s\n", dirpath);
			exit( - 1);
		}
	}
	return ;
}

/* Parse the command line arguments */
void parseArguments(int argc, char * * argv, char * * input, char * * events, char * * ipts, char * devicePref, short * verifyResults)
{
	{
		int i = 2;
		#pragma loop name parseArguments#0 
		for (; i<argc; i ++ )
		{
			if (strcmp(argv[i], "-d")==0)
			{
				/* Event dump found */
				if (i==(argc-1))
				{
					printf("Usage: -e Needs directory path\n");
					exit( - 1);
				}
				devicePref[0]=argv[i+1][0];
				i ++ ;
				continue;
			}
			if (strcmp(argv[i], "-e")==0)
			{
				/* Event dump found */
				if (i==(argc-1))
				{
					printf("Usage: -e Needs directory path\n");
					exit( - 1);
				}
				( * events)=argv[i+1];
				i ++ ;
				continue;
			}
			if (strcmp(argv[i], "-i")==0)
			{
				/* Input found */
				if (i==(argc-1))
				{
					printf("Usage: -i Needs directory path\n");
					exit( - 1);
				}
				( * input)=argv[i+1];
				i ++ ;
				continue;
			}
			if (strcmp(argv[i], "-l")==0)
			{
				/* Ipts dump found */
				if (i==(argc-1))
				{
					printf("Usage: -l Needs directory path\n");
					exit( - 1);
				}
				( * ipts)=argv[i+1];
				i ++ ;
				continue;
			}
			if (strcmp(argv[i], "-n")==0)
			{
				/* Don't use OpenCL images */
				setUsingImages(0);
				continue;
			}
			if (strcmp(argv[i], "-v")==0)
			{
				/* Verify results */
				( * verifyResults)=1;
				continue;
			}
		}
	}
	return ;
}

/* This function that takes a positive integer 'value' and returns */
/* the nearest multiple of 'multiple' (used for padding columns) */
unsigned int roundUp(unsigned int value, unsigned int multiple)
{
	unsigned int remainder = value%multiple;
	/* Make the value a multiple of multiple */
	if (remainder!=0)
	{
		value+=(multiple-remainder);
	}
	return value;
}

/* Concatenate two strings and return a pointer to the new string */
char *smartStrcat(char * str1, char * str2)
{
	char * newStr = (void * )0;
	newStr=((char * )alloc(((strlen(str1)+strlen(str2))+1)*sizeof (char)));
	strcpy(newStr, str1);
	strcat(newStr, str2);
	return newStr;
}

/* Set the value of using images to true if they are being */
/* used, or false if they are not */
void setUsingImages(short val)
{
	usingImages=val;
	return ;
}

/* Return whether or not images are being used */
short isUsingImages()
{
	return usingImages;
}
