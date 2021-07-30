include ../../common/make.config

#C compiler
# CC = g++
CC = gcc

# SRC = clutils.cpp gaussianElim.cpp utils.cpp

CC_FLAGS = -O3
LDLAGS  = -lm -fPIC

# CC_FLAGS = -g -O3

# EXE = gaussian

gaussian_host: gaussianElim.o clutils.o utils.o
	$(CC) $(KERNEL_DIM) $(CC_FLAGS) --std=c99  gaussianElim.c clutils.o utils.o -o gaussian_host $(LDLAGS) -g `dpu-pkg-config --cflags --libs dpu`

%.o: %.[ch]
	$(CC) $(CC_FLAGS) $< -c

# release:
# 	$(CC) $(KERNEL_DIM) $(CC_FLAGS) -lOpenCL  -fopenmp $(SRC) -o $(EXE) -I$(OPENCL_INC) -L$(OPENCL_LIB)


clutils.o: clutils.c
	$(CC) $(CC_FLAGS) clutils.c -c $(LDLAGS)

utils.o: utils.c
	$(CC) $(CC_FLAGS) utils.c -c $(LDLAGS)

clean: 
	rm -f $(EXE) 
