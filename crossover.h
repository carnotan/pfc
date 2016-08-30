#ifndef CROSSOVER_H
#define CROSSOVER_H

float get_beta(float u, float eps, float d_index);
int sbx(Solution *parent1, Solution *parent2, float
        cross_prob, float d_index, float *upper_b, float *lower_b, float epsilon);
int recombine(thrust::host_vector <Solution> *mating_pool, size_t pop_size, float
        cross_prob, float d_index, float *upper_b, float *lower_b, float epsilon);

#endif