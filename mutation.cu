#include "data_structures.h"
#include <thrust/host_vector.h>

/**
 * Función auxiliar que recalcula os parámetros da mutación (índice e taxa de 
 * mutación).
 * 
 * @param index Punteiro ao índice da función de distribución de probabilidade
 *  polinomial.
 * @param rate Punteiro á taxa de mutación.
 * @param g_size Número de variables do cromosoma.
 * @param ex_count Ciclos do algoritmo executados.
 * @param max_exec Máximo de ciclos que se poden executar
 * @param d_index Índice da función de distribución de probabilidade polinomial 
 * orixinal.
 */
void recalculate_parameters(float *index, float *rate, float g_size, 
        int ex_count, int max_exec, float d_index) {

    *rate = 1.0f / g_size + ((float) ex_count / (float) max_exec)*
            (1.0f - (1.0f / g_size));
    *index = d_index + ex_count;
    
}

/**
 * Función que obtén o delta para a mutación polinomial.
 * 
 * @param u  Parámetro de distribución uniforme u[0,1]
 * @param delta_l Parámetro delta_lower.
 * @param delta_u Parámetro delta_upper.
 * @param d_index Punteiro ao índice da función de distribución de probabilidade
 *  polinomial.
 * @return O valor delta para a mutación .
 */
float get_delta(float u, float delta_l, float delta_u, float *d_index) {
    float delta, aa;

    if (u >= 1.0f - 1.0e-7) delta = delta_u;
    else if (u <= 0.0f + 1.0e-7) delta = delta_l;
    else {
        if (u <= 0.5f) {
            aa = 2.0f * u + (1.0f - 2.0f * u) * powf((1.0f + delta_l), 
                    (*d_index + 1.0f));
            delta = powf(aa, (1.0f / (*d_index + 1.0f))) - 1.0f;
        } else {
            aa = 2.0f * (1.0f - u) + 2.0f * (u - 0.5f) * powf((1.0f - delta_u),
                    (*d_index + 1.0f));
            delta = 1.0f - powf(aa, (1.0 / (*d_index + 1.0f)));
        }
    }
    if (delta < -1.0f || delta > 1.0f) {
        printf("Erro na mutación!! delta = %lf\n", delta);
        exit(-1);
    }
    return (delta);
}

/**
 * Función que realiza a mutación de un individuo usando unha función polinomial
 * para perturbar os valores dos xenes nas proximidades do pai (mutación paterno-
 * céntrica) . Implementación con límites inferior e superior proposta por K.Deb
 * e R.Agrawal.
 * 
 * @param parent Punteiro á estrutura solution a mutar.
 * @param d_index Punteiro ao índice da función de distribución de probabilidade polinomial. 
 * @param upper_b Punteiro ao array de límites superiores.
 * @param lower_b Punteiro ao array de límites inferiores.
 * @param mut_rate Punteiro á taxa de mutación.
 */
void real_mutation(Solution* parent, float *d_index, float *upper_b,
        float *lower_b, float *mut_rate) {

    float delta, u;
    float delta_l, delta_u;
    int i;

    for (i = 0; i < 4; i++) {
        if (((float) rand() / (float) RAND_MAX) <= *mut_rate) {
            delta_l = (lower_b[i] - parent->chromosome[i]) / (upper_b[i] - 
                    lower_b[i]);
            if (delta_l < -1.0f)
                delta_l = -1.0f;
            delta_u = (upper_b[i] - parent->chromosome[i]) / (upper_b[i] - 
                    lower_b[i]);
            if (delta_u > 1.0f)
                delta_u = 1.0f;
            if (-1.0f * delta_l < delta_u)
                delta_u = -1.0f * delta_l;
            else delta_l = -1.0f * delta_u;
            u = ((float) rand() / (float) RAND_MAX);
            delta = get_delta(u, delta_l, delta_u, d_index)
                    * (upper_b[i] - lower_b[i]);
            parent->chromosome[i] = parent->chromosome[i] + delta;
            if (parent->chromosome[i] < lower_b[i])
                parent->chromosome[i] = lower_b[i];
            if (parent->chromosome[i] > upper_b[i])
                parent->chromosome[i] = upper_b[i];
        }
    }

}

/**
 * Función auxiliar que recorre toda a poboación realizando a mutación individuo 
 * a individuo.
 * 
 * @param pop_size Tamaño da poboación.
 * @param mating_pool Punteiro a host_vector<Solution> que contén a piscina de 
 * reprodución.
 * @param mut_rate Punteiro á taxa de mutación.
 * @param d_index Punteiro ao índice da función de distribución de probabilidade
 *  polinomial. 
 * @param upper_b Punteiro ao array de límites superiores.
 * @param lower_b Punteiro ao array de límites inferiores.
 */
void mutation(size_t pop_size, thrust::host_vector <Solution> *mating_pool,
        float * mut_rate, float *d_index, float *upper_b, float *lower_b) {

    for (int i = 0; i < pop_size; i++) {
        real_mutation(&mating_pool->operator[](i), d_index, upper_b, lower_b, 
                mut_rate);
    }
    
}
