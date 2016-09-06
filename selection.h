#ifndef SELECTION_H
#define SELECTION_H

int get_winner(std::vector <Solution> pop, size_t tour_size, 
        int * players);
int tournament_selection(std::vector <Solution> tour_pop, 
        std::vector <Solution>*mating_pool, size_t pop_size, 
        size_t tour_size);

#endif