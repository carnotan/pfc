#include "data_structures.h"
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "auxiliares.h"

/**
 * Función que determina o gañador dos torneos de selección para a piscina de
 * reprodución.
 * 
 * @param pop host_vector<Solution> que contén a poboación.
 * @param tour_size Tamaño do torneo.
 * @param players Punteiro ao array que contén o índice (dentro da poboación) 
 * dos participantes no torneo.
 * @return O índice do gañador.
 */
int get_winner(std::vector <Solution> pop, size_t tour_size, 
        int * players) {

    int index = players[0];
    float a = pop[players[0]].fitness;

    for (int i = 1; i < tour_size; i++) {
        if (pop[players[i]].fitness > a) {
            a = pop[players[i]].fitness;
            index = players[i];
        }
    }
    return index;

}

/**
 * Función que xera a piscina de reproducción a partir da poboación realizando 
 * torneos sucesivos de tamaño tour_size. Selección determinística por fitness.
 * 
 * @param tour_pop host_vector<Solution> que contén a poboación.
 * @param mating_pool Punteiro a host_vector<Solution> que contén a piscina de 
 * reproducción.
 * @param pop_size Tamaño da poboación.
 * @param tour_size Tamaño do torneo.
 * @return 0 se a función remata correctamenet, código do erro noutro caso.
 */
int tournament_selection(std::vector <Solution> tour_pop,
        std::vector <Solution>*mating_pool, size_t pop_size,
        size_t tour_size) {

    int winner;
    int result;
    int *aux = (int *) malloc(tour_size * sizeof (int));
    
    if (aux==NULL){
        printf("Non hai memoria suficiente para facer malloc\n");
        printf("Erro realizando a selección\n");
        return (-1);
    }
    for (int i = 0; i <= pop_size - tour_size; i++) { 
        if ((result = shuffle(aux, pop_size - i, tour_size)) == -1) {
            return result;
        }
        winner = get_winner(tour_pop, tour_size, aux);
        mating_pool->operator[](i) = tour_pop[winner];
        tour_pop[winner] = tour_pop[pop_size - i - 1];
        tour_pop[pop_size - i - 1] = mating_pool->operator[](i);
    }
    //Os derradeiros elementos son copiados directamente á piscina.
    for (int i = pop_size - tour_size + 1; i < pop_size; i++) {
        mating_pool->operator[](i) = tour_pop[i - pop_size + tour_size - 1];
    }
    return 0;
    
}
