#include "data_structures.h"
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/**
 * Función para o cálculo do parámetro beta necesario para o SBX.
 * 
 * @param u Parámetro de distribución uniforme u[0,1]
 * @param eps Épsilon da máquina.
 * @param d_index Índice da distribución de probabilidade.
 * @return O valor beta para o SBX.
 */
float get_beta(float u, float eps, float d_index) {

    float beta;

    if (1.0f - u < eps)
        u = 1.0f - eps;
    if (u < 0.0f)
        u = 0.0f;
    if (u < 0.5f)
        beta = powf(2.0f * u, (1.0 / (d_index + 1.0f)));
    else
        beta = powf((0.5f / (1.0f - u)), (1.0f / (d_index + 1.0f)));
    return beta;

}

/**
 * Simulated Binary Crossover proposta por K.Deb e R.Agrawal (1995). Variante 
 * con límite superior e inferior. Realiza o cruce de dous individuos empregando
 * unha distribución de probabilidade para simular cruces de 1 punto en 
 * variables binarias. O cruce realízase con probabilidade xeral cross_prob e 
 * con probabilidade de cruce 0.5 .
 * 
 * @param parent1 Estructura solution do pai 1.
 * @param parent2 Estructura solution do pai 2.
 * @param cross_prob Probabilidade de cruce (0.5).
 * @param d_index Índice para a función de probabilidade polinomial.
 * @param upper_b Límite superior para o cruce.
 * @param lower_b Límite inferior para o cruce.
 * @param epsilon Épsilon de máquina.
 * @result 0 se se fai correctamente, -1 noutro caso.
 */
int sbx(Solution *parent1, Solution *parent2, float
        cross_prob, float d_index, float *upper_b, float *lower_b,
        float epsilon) {

    float beta;
    float u;
    float distance, difference;
    float x_mean;
    float umax;
    float alpha;
    int i;
    Solution *temp = (Solution*) malloc(sizeof (Solution));

    if (temp == NULL)
        return (-1);

    if (((float) rand() / (float) RAND_MAX) < cross_prob) {
        for (i = 0; i < 4; i++) {
            if (((float) rand() / (float) RAND_MAX) <= 0.5f) {
                if (parent2->chromosome[i] < parent1->chromosome[i]) {
                    memcpy(temp, parent2, sizeof (Solution));
                    memcpy(parent2, parent1, sizeof (Solution));
                    memcpy(parent1, temp, sizeof (Solution));
                }
                x_mean = (parent1->chromosome[i] + parent2->chromosome[i])*0.5f;
                difference = parent2->chromosome[i] - parent1->chromosome[i];
                if (difference > epsilon) {
                    if ((parent1->chromosome[i] - lower_b[i])< (upper_b[i] -
                            parent2->chromosome[i]))
                        distance = parent1->chromosome[i] - lower_b[i];
                    else
                        distance = upper_b[i] - parent2->chromosome[i];
                    if (distance < 0.0f) distance = 0.0f;
                    alpha = 1.0f + (2.0f * distance / difference);
                    umax = 1.0f - (0.5f / powf(alpha, (d_index + 1.0f)));
                    u = umax * ((float) rand() / (float) RAND_MAX);
                } else
                    u = ((float) rand() / (float) RAND_MAX);
                beta = get_beta(u, epsilon, d_index);
                parent1->chromosome[i] = x_mean - beta * 0.5f * difference;
                parent2->chromosome[i] = x_mean + beta * 0.5f * difference;
            }
        }
    }
    return 0;

}

/**
 * Función auxiliar. Recorre a poboación e realiza o cruce 2 a 2 , depositándoos 
 * na poboación de saída.
 * 
 * @param mating_pool Poboación de saída / piscina de reprodución.
 * @param pop_size Tamaño da poboación.
 * @param cross_prob Probabilidade de cruce.
 * @param d_index Índice para a función de probabilidade polinomial.
 * @param upper_b Límite superior para o cruce.
 * @param lower_b Límite inferior para o cruce.
 * @param epsilon Épsilon de máquina.
 */
int recombine(std::vector <Solution> *mating_pool, size_t pop_size, float
        cross_prob, float d_index, float *upper_b, float *lower_b,
        float epsilon) {

    int result;

    for (int i = 0; i < pop_size; i += 2) {
        if ((result = sbx(&mating_pool->operator[](i), 
                &mating_pool->operator[](i + 1), cross_prob, d_index, upper_b,
                        lower_b, epsilon)) == -1)
            return result;
    }
    return result;
    
}
