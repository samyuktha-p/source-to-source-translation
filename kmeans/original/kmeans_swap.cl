#ifndef FLT_MAX
#define FLT_MAX 3.40282347e+38
#endif

void kmeans_swap(float  *feature, float  *feature_swap, int npoints, int nfeatures){

	unsigned int tid = get_global_id(0);
	//for(int i = 0; i <  nfeatures; i++)
	//	feature_swap[i * npoints + tid] = feature[tid * nfeatures + i];
    //Lingjie Zhang modificated at 11/05/2015
    if (tid < npoints){
	    for(int i = 0; i <  nfeatures; i++)
		    feature_swap[i * npoints + tid] = feature[tid * nfeatures + i];
    }
    // end of Lingjie Zhang's modification
} 
