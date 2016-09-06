#ifndef AUXILIARES_H
#define AUXILIARES_H


int reseed(void);
int shuffle(int *array, size_t N, size_t M);

int is_converged(std::vector<Solution> population, size_t pop_size,
        float *previous, int * not_improved, bool fast_convergence, 
        float min_growth);
void write_solution(Solution s) ;

#endif
