include ../../common/make.config

#C compiler
# CC = g++
CC = gcc

# SRC = clutils.cpp gaussianElim.cpp utils.cpp

CC_FLAGS = -O3
LDLAGS  = -lm -fPIC

# CC_FLAGS = -g -O3

# EXE = gaussian
kmeans_host: cluster.o getopt.o read_input.o kmeans_clustering.o rmse.o
	$(CC) $(KERNEL_DIM) $(CC_FLAGS) --std=c99  kmeans.c cluster.o getopt.o read_input.o kmeans_clustering.o rmse.o -o kmeans_host $(LDLAGS) -g `dpu-pkg-config --cflags --libs dpu`

%.o: %.[ch]
	$(CC) $(CC_FLAGS) $< -c

# release:
# 	$(CC) $(KERNEL_DIM) $(CC_FLAGS) -lOpenCL  -fopenmp $(SRC) -o $(EXE) -I$(OPENCL_INC) -L$(OPENCL_LIB)


%.o: %.[ch]
	$(CC) $(CC_FLAGS) $< -c

cluster.o: cluster.c 
	$(CC) $(CC_FLAGS) cluster.c -c $(LDLAGS)

getopt.o: getopt.c 
	$(CC) $(CC_FLAGS) getopt.c -c $(LDLAGS)

kmeans.o: kmeans.c 
	$(CC) $(CC_FLAGS) read_input.c -c $(LDLAGS)

rmse.o: rmse.c
	$(CC) $(CC_FLAGS) rmse.c -c $(LDLAGS)

clean:
	rm -f *.o *~ kmeans *.linkinfo

#clean: 
#	rm -f $(EXE) 
