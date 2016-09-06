#ifndef DEBUG_H
#define DEBUG_H

void write_eliminated(std::vector<Point> cloud, size_t cloud_size, const char *str, std::vector<int> vector/*,std::ios::openmode mode*/);
void write_eliminated_serie(std::vector<Point> cloud, size_t eliminados, const char *str /*,std::ios::openmode mode*/);
static void show_chromosome(Solution * s);
static void show_population(std::vector<Solution>genome, size_t pop_size);
static void show_cloud(std::vector<Point> cloud, size_t cloud_size);
 void write_cloud(std::vector<Point> cloud, size_t cloud_size, const char *str/*,std::ios::openmode mode*/);
void write_population(std::vector<Solution> population, size_t pop_size, const char *str/*,std::ios::openmode mode*/);
void write_vector (std::vector<int>,size_t size,const char *str);
void show_vector(std::vector <int> vector, size_t size);

#endif