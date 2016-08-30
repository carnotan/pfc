#include "data_structures.h"
#include "thrust/host_vector.h"
#include <float.h>

/**
 * Función que obtén o fitness asociado a un individuo da pobación. Obtén a 
 * distancia euclídea de un punto respecto do plano hipotetizado na solución. 
 * 
 * @param t Umbral do plano (threshold).
 * @param r Rexión do plano (region).
 * @param eps Épsilon de máquina.
 * @param cloud host_vector<Solution> que contén a nube de puntos.
 * @param b Base do logaritmo da función de fitness.
 * @param p Punteiro a estrutura Solution cos parámetros do plano a evaluar.
 * @param cloud_size Tamaño da nube.
 * @return O valor de fitness para este punto.
 */
float fitness(float t, float r, float eps, thrust::host_vector<Point> cloud,
        float b, Solution * p, size_t cloud_size) {

    float d = 0.0f;
    float k = 0.0f;
    int i;
    int l = 0;
    float fitness = 0.0f;
    float max_log;
    p->points_in_region = 0;

    b = powf((r / t), (1.0f / r));
    max_log = (log10f(r) / log10f(b))-(log10f(t) / log10f(b));
    for (i = 0; i < cloud_size; i++) {
        d = fabsf(cloud[i].coordinates[0] * p->chromosome[0] +
                cloud[i].coordinates[1] * p->chromosome[1] +
                cloud[i].coordinates[2] * p->chromosome[2] + p->chromosome[3]);
        k = sqrtf(powf(p->chromosome[0], 2.0f) + powf(p->chromosome[1], 2.0f)
                + powf(p->chromosome[2], 2.0f));
        if (k <= 5 * FLT_EPSILON)
            d = d / eps;
        else
            d = d / k;
        if (d > 5 * FLT_EPSILON && d <= r) {
            p->points_in_region = p->points_in_region + 1;
            fitness = fitness + (log10f(d) / log10f(b))-(log10f(t) / log10f(b))
                    - max_log;
        }
        if (d > 5 * FLT_EPSILON && d <= t)
            l++;
        if (d < 5 * FLT_EPSILON)
            p->points_in_region = p->points_in_region + 1;
    }
    if (p->points_in_region == 0)
        fitness = -INFINITY;
    else
        fitness = -fitness / p->points_in_region * l;
    p->points_fitted = l;
    return fitness;

}

/**
 * Función auxiliar que chama á función de fitness para reevaluar a calidade da
 * poboación.
 *
 * @param t Umbral do plano (threshold).
 * @param r Rexión do plano (region).
 * @param eps Epsilon de máquina.
 * @param cloud host_vector<Point> que contén a nube de puntos.
 * @param b Base do logaritmo da función de calidade. 
 * @param population Punteiro a host_vector<Solution> que contén a poboación.
 * @param cloud_size Tamaño da nube de puntos.
 * @param pop_size Tamaño da poboación.
  */
void evaluate_population(float t, float r, float eps,
        thrust::host_vector<Point> cloud, float b,
        thrust::host_vector<Solution> *population, size_t cloud_size,
        size_t pop_size) {

    for (int i = 0; i < pop_size; i++) {
        population->operator[](i).fitness = fitness(t, r, eps, cloud, b,
                &population->operator[](i), cloud_size);
    }

}