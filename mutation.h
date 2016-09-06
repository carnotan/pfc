#ifndef MUTATION_H
#define MUTATION_H

float get_delta(float u, float delta_l, float delta_u, float *d_index);
void real_mutation(Solution* parent, float *d_index, float *upper_b, 
        float *lower_b, float *mut_rate);
void mutation(size_t pop_size, std::vector <Solution> *mating_pool, 
        float * mut_rate, float *d_index, float *upper_b, float *lower_b);
void recalculate_parameters(float *index, float *rate, float g_size, 
        int ex_count, int max_exec, float d_index);

#endif
