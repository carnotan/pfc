#ifndef SELECTION_H
#define SELECTION_H

int get_winner(thrust::host_vector <Solution> pop, size_t tour_size, 
        int * players);
int tournament_selection(thrust::host_vector <Solution> tour_pop, 
        thrust::host_vector <Solution>*mating_pool, size_t pop_size, 
        size_t tour_size);

#endif