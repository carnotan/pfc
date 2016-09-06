
#include "data_structures.h"
#include <cstddef>
#include <vector>
#include <float.h>
#include <math.h>

/**
 * Función que elimina da nube de puntos aqueles que pertencen ao plano ao que 
 * converxeu o algoritmo.
 * 
 * @param p Estructura Solution que contén os parámetros do plano atopado.
 * @param cloud_size Tamaño da nube de puntos.
 * @param cloud Punteiro a std::vector<Point> que contén a nube de puntos.
 * @param t Umbral do plano.
 * @param eps Épsilon de máquina.
 */
void eliminate(Solution p, size_t *cloud_size, std::vector<Point> *cloud,
        float t, float eps) {

    int deleted = 0;
    float d, k;

    for (int i = 0; i < *cloud_size && deleted < p.points_fitted; i++) {
        d = fabsf(cloud->operator[](i).coordinates[0] * p.chromosome[0] + 
        cloud->operator[](i).coordinates[1] * p.chromosome[1] + 
        cloud->operator[](i).coordinates[2] * p.chromosome[2] + p.chromosome[3]);
        k = sqrtf(powf(p.chromosome[0], 2.0f) + powf(p.chromosome[1], 2.0f) 
                + powf(p.chromosome[2], 2.0f));
        if (k >= 0.0f && k <= 5 * FLT_EPSILON)
            d = d / eps;
        else
            d = d / k;
        if (d > 5 * FLT_EPSILON && d <= t) {
            deleted++;
            std::iter_swap(cloud->begin() + i, cloud->end() - deleted);
        }
    }
    cloud->erase(cloud->end() - deleted, cloud->end());
    *cloud_size = *cloud_size - deleted;
    
}