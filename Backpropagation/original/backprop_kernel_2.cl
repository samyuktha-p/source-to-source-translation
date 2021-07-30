#define THREADS 256
#define WIDTH 16  
#define HEIGHT 16 
#define ETA 0.3f       
#define MOMENTUM 0.3f  

#ifndef _BACKPROP_CUDA_KERNEL_H_2
#define _BACKPROP_CUDA_KERNEL_H_2
#define WM(i, j)   weight_matrix[(j) + (i) * WIDTH]




void  bpnn_adjust_weights_ocl( float * delta,   
										 int hid,         
										float * ly,      
										 int in,          
										float * w,       
										float * oldw)  									
{
   
   int by = get_group_id(1);
   int tx = get_local_id(0);
   int ty = get_local_id(1);
	
   // int index =  ( hid + 1 ) * HEIGHT * by + ( hid + 1 ) * ty + tx + 1 + ( hid + 1 ) ;  
   // int index_y = HEIGHT * by + ty + 1;
   // int index_x = tx + 1;

   w[( hid + 1 ) * HEIGHT * by + ( hid + 1 ) * ty + tx + 1 + ( hid + 1 )] += ((ETA * delta[tx + 1] * ly[HEIGHT * by + ty + 1]) + (MOMENTUM * oldw[( hid + 1 ) * HEIGHT * by + ( hid + 1 ) * ty + tx + 1 + ( hid + 1 )]));
   oldw[( hid + 1 ) * HEIGHT * by + ( hid + 1 ) * ty + tx + 1 + ( hid + 1 )] = ((ETA * delta[tx + 1] * ly[HEIGHT * by + ty + 1]) + (MOMENTUM * oldw[( hid + 1 ) * HEIGHT * by + ( hid + 1 ) * ty + tx + 1 + ( hid + 1 )]));

   barrier(CLK_LOCAL_MEM_FENCE);

   if (ty == 0 && by ==0){
	w[tx + 1] += ((ETA * delta[tx + 1]) + (MOMENTUM * oldw[tx + 1]));
	oldw[tx + 1] = ((ETA * delta[tx + 1]) + (MOMENTUM * oldw[tx + 1]));
   }

}
#endif 
