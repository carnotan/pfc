#ifndef ELIMINATE_H
#define ELIMINATE_H 

void eliminate(Solution p, size_t *cloud_size, thrust::host_vector<Point> *cloud,
        float t, float eps);

#endif

