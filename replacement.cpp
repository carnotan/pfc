#include "data_structures.h"
#include <vector>
#include <algorithm>

/**
 * Función que implementa o algoritmo de reemplazo. Óptase por unha estratexia 
 * de elitismo puro para garantizar a converxencia da poboación . Dadas a 
 * poboación orixinal e a piscina de reprodución (mutada e cruzada) de xeración
 *  N, a xeración N+1 obtense mezclando ambas poboacións e quedándonos coa 
 * metade mellor (en función do fitness) de ambas. 
 * 
 * @param pop Punteiro a host_vector<Solution> que contén a poboación de 
 * xeración N
 * @param mating_pool Punteiro a host_vector<Solution> que contén a piscina de 
 * reprodución de xeración N .
 * @param pop_size Tamaño da poboación.
 */
void replacement(std::vector<Solution> *pop, 
        std::vector<Solution> mating_pool, size_t pop_size) {
    
    std::vector <Solution> aux(2 * pop_size);
    std::copy(pop->begin(), pop->end(), aux.begin());
    std::copy(mating_pool.begin(), mating_pool.end(), aux.begin() + pop_size);
    std::sort(aux.begin(), aux.end(), SOLUTIONCmp());
    std::copy(aux.begin(), aux.begin() + pop_size, pop->begin());
    
}