#ifndef AUXILIARES_H
#define AUXILIARES_H

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

int reseed(void);
int shuffle(int *array, size_t N, size_t M);

int is_converged(thrust::host_vector <Solution> population, size_t pop_size,
        float *previous, int * not_improved, bool fast_convergence,
        float min_growth);
void write_solution(Solution s);
float average_fitness(thrust::host_vector<Solution> population,
        size_t pop_size) ;
inline void gpuAssert(cudaError_t code, const char *file, int line,
        bool abort = true) {
    if (code != cudaSuccess) {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file,
                line);
        if (abort) exit(code);
    }
}
#endif
