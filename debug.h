#ifndef DEBUG_H
#define DEBUG_H

void write_eliminated(thrust::host_vector<Point> cloud, size_t cloud_size, const char *str, thrust::device_vector<int> vector/*,std::ios::openmode mode*/);
void write_eliminated_serie(thrust::host_vector<Point> cloud, size_t eliminados, const char *str /*,std::ios::openmode mode*/);
static void show_chromosome(Solution * s);
static void show_population(thrust::host_vector<Solution>genome, size_t pop_size);
static void show_cloud(thrust::host_vector<Point> cloud, size_t cloud_size);
 void write_cloud(thrust::host_vector<Point> cloud, size_t cloud_size, const char *str/*,std::ios::openmode mode*/);
void write_population(thrust::host_vector<Solution> population, size_t pop_size, const char *str/*,std::ios::openmode mode*/);
void write_vector (thrust::device_vector<int>,size_t size,const char *str);
void show_vector(thrust::device_vector <int> vector, size_t size);
int contar(thrust::device_vector<int> vector);

#endif