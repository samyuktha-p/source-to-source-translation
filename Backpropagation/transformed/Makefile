#include ../../common/make.config

# C compiler
CC = gcc
#CC_FLAGS = -g -O3 -lm -ldl -Wall -Wpointer-arith -finline-functions -ffast-math -funroll-all-loops

CC_FLAGS = -O3
LDLAGS  = -lm -fPIC
#backprop: backprop.o imagenet.o facetrain.o
#	$(CC) $(CC_FLAGS) backprop_ocl.cpp backprop.o imagenet.o facetrain.o -o backprop -I$(OPENCL_INC) -L$(OPENCL_LIB) -lOpenCL
backprop_host: backprop.o imagenet.o facetrain.o
	$(CC) $(CC_FLAGS) --std=c99  backprop_ocl.c backprop.o imagenet.o facetrain.o -o backprop_host $(LDLAGS) -g `dpu-pkg-config --cflags --libs dpu`

%.o: %.[ch]
	$(CC) $(CC_FLAGS) $< -c
	
backprop.o: backprop.c 
	$(CC) $(CC_FLAGS) backprop.c -c $(LDLAGS)
	
facetrain.o: facetrain.c
	$(CC) $(CC_FLAGS) facetrain.c -c $(LDLAGS)

imagenet.o: imagenet.c
	$(CC) $(CC_FLAGS) imagenet.c -c $(LDLAGS)

clean:
	rm -f *.o *~ backprop *.linkinfo
