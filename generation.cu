#include "data_structures.h"
#include "thrust/host_vector.h"
#include "auxiliares.h"

/**
 * Función que obtén os parámetros da ecuación xeral dun plano que pasa por 3
 *  puntos.
 * 
 * @param p1 Punto 1 da nube.
 * @param p2 Punto 2 da nube.
 * @param p3 Punto 3 da nube.
 * @param s Punteiro á estrutura Solution de saída cos parámetros do plano 
 * xerado. 
 */
void get_solution(Point p1, Point p2, Point p3, Solution *s) {

    s->chromosome[0] = (((p2.coordinates[1])-(p1.coordinates[1]))*
            ((p3.coordinates[2])-(p1.coordinates[2])))-(((p3.coordinates[1])-
            (p1.coordinates[1]))*((p2.coordinates[2])-(p1.coordinates[2])));
    s->chromosome[1] = ((-1.0f)*((p2.coordinates[0])-(p1.coordinates[0]))*
            ((p3.coordinates[2])-(p1.coordinates[2])))+(((p3.coordinates[0])-
            (p1.coordinates[0]))*((p2.coordinates[2])-(p1.coordinates[2])));
    s->chromosome[2] = (((p2.coordinates[0])-(p1.coordinates[0]))*
            ((p3.coordinates[1])-(p1.coordinates[1])))-(((p3.coordinates[0])-
            (p1.coordinates[0]))*((p2.coordinates[1])-(p1.coordinates[1])));
    s->chromosome[3] = ((-1.0f) * s->chromosome[0] * p1.coordinates[0] -
            s->chromosome[1] * p1.coordinates[1] - s->chromosome[2] 
            * p1.coordinates[2]);

}

/**
 * Función que normaliza os valores dos cromomas en vectores de lonxitude 1 de
 * acordo á formula : v=(x, y, z) 
 * vu = v / |v| = ( x /|v|,  y/|v|,  z/|v| )
 * 
 * @param population Punteiro ao host_vector<Solution> que contén á poboación.
 * @param pop_size Tamaño da poboación.
 */

void normalize(thrust::host_vector<Solution> *population, size_t pop_size) {
    float norm;

    for (int i = 0; i < pop_size; i++) {
        norm = sqrtf(powf(population->operator[](i).chromosome[0], 2.0f) + 
                powf(population->operator[](i).chromosome[1], 2.0f) + 
                powf(population->operator[](i).chromosome[2], 2.0f) + 
                powf(population->operator[](i).chromosome[3], 2));
        for (int j = 0; j < 4; j++) {
            population->operator[](i).chromosome[j] = 
                    population->operator[](i).chromosome[j] / norm;
        }
    }

}

/**
 * Función que crea unha poboación de tamaño N tomando n*3 puntos non repetidos 
 * da nube.
 * 
 * @param pop_size Tamaño da poboación.
 * @param cloud_size Tamaño da nube de puntos.
 * @param population Punteiro a un host_vector<Solution> que contén a poboación.
 * @param cloud Punteiro a un host_vector<Point> que contén a nube de puntos.
 * @return 0 se a poboación se creou correctamente, o código de erro noutro
 *  caso.
 */
int generate_population(size_t pop_size, size_t cloud_size, 
        thrust::host_vector<Solution> *population, 
        thrust::host_vector<Point> *cloud) {

    int i,result;
    
    if (cloud_size < 3 * pop_size) {
        printf("Tamaño insuficiente da nube para xenerar a poboación\n");
        return (1);
    }
    int *array = (int*) malloc(3 * pop_size * sizeof (int));
    if (array == NULL) {
        printf("Non hai memoria suficiente para facer malloc\n");
        printf("Erro xerando a poboación\n");
        return (-1);
    }
    if ((result = shuffle(array, cloud_size, pop_size * 3)) == -1) {
        return (result);
    }
    for (i = 0; i < pop_size; i++) {
        get_solution(cloud->operator[](array[3 * i]),
                cloud->operator[](array[3 * i + 1]), 
                        cloud->operator[](array[3 * i + 2]),
                                &population->operator[](i));
    }
    return 0;
    
}
