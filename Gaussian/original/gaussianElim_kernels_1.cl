typedef struct latLong
    {
        float lat;
        float lng;
    } LatLong;

void Fan1(float *m_dev,
                  float *a_dev,
                  float *b_dev,
                  long size,
                  long t) {
    int globalId = get_global_id(0);
                              
    if (globalId < size-1-t) {
         //*(m_dev + size * (globalId + t + 1)+t) = *(a_dev + size * (globalId + t + 1) + t) / *(a_dev + size * t + t);    
			m_dev[size * (globalId + t + 1)+t] = a_dev[size * (globalId + t + 1) + t] / a_dev [size * t + t];        
    }
}
