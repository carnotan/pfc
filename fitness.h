#ifndef FITNESS_H
#define FITNESS_H

__global__ void
__launch_bounds__(128, 8)
fitness(float t, float r, float eps, Point *cloud, Solution *s, int * fitted,
        int * region, float * fit, float b, float max_log);

__device__ int getGlobalIdx_1D_1D();
int getBlockSize(int original, int max);
void evaluate_population_cuda(float t, float r, float eps, Point *p_d_cloud,
        thrust::host_vector<Solution> *population, size_t cloud_size,
        size_t pop_size, thrust::device_vector <int> *fitted, 
        thrust::device_vector <int> *region, thrust::device_vector<float>*fit);

#endif