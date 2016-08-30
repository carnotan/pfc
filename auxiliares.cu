#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <thrust/host_vector.h>
#include "data_structures.h"
#include <fstream>

//Estructuras e variables empregadas na xeración de números aleatorios

static struct {
    int which;
    time_t t;
    clock_t c;
    int counter;
} entropy = {0, (time_t) 0, (clock_t) 0, 0};

static unsigned char * pSeed = (unsigned char *) (&entropy + 1);
static int accSeed = 0;

/**
 * Función auxiliar que devolve unha nova semente que pasar as funcións srand
 * ou srand48. Permite ser chamada múltiples veces por segundo.
 * 
 * @return Un enteiro (semente).
 */
int reseed(void) {
    if (pSeed == ((unsigned char *) (&entropy + 1))) {
        switch (entropy.which) {
            case 0:
                entropy.t += time(NULL);
                accSeed ^= entropy.t;
                break;
            case 1:
                entropy.c += clock();
                break;
            case 2:
                entropy.counter++;
                break;
        }
        entropy.which = (entropy.which + 1) % 3;
        pSeed = (unsigned char *) &entropy.t;
    }
    accSeed = ((accSeed * (UCHAR_MAX + 2U)) | 1) + (int) *pSeed;
    pSeed++;
    return accSeed;
}

/**
 * Función auxiliar que devolve un array de M enteiros non repetidos
 * seleccionados aleatoriamente de entre N elementos equiprobables. Segue a 
 * implementación do algoritmo de Floyd-Warshall
 * https://en.wikipedia.org/wiki/Floyd%E2%80%93Warshall_algorithm
 * 
 * @param array Array de enteiros devoltos.
 * @param N tamaño da poboación.
 * @param M tamaño da mostra tomada.
 * @return 0 se se xenera o array correctamente, erro noutro caso.
 */
int shuffle(int *array, size_t N, size_t M) {
    unsigned char *is_used = (unsigned char *) malloc(N * sizeof (unsigned char));
    int in, im;

    if (is_used == NULL) {
        printf("Non hai memoria suficiente para facer malloc\n");
        return (-1);
    }
    memset(is_used, 0x0, N * sizeof (unsigned char));
    im = 0;
    srand(reseed());
    for (in = N - M; in < N && im < M; ++in) {
        int r = rand() % (in + 1);
        if (is_used[r])
            r = in;
        assert(!is_used[r]);
        array[im++] = r;
        is_used[r] = 1;
    }
    assert(im == M);
    return (0);
}

/**
 * Calcula se a poboación converxeu dacordo ao criterio do algoritmo .
 * 
 * @param population host_vector de estructuras Solution.
 * @param pop_size tamaño da poboación
 * @return 1 se a poboación converxeu, 0 noutro caso
 */
int is_converged(thrust::host_vector<Solution> population, size_t pop_size) {

    int a = population[0].points_fitted;

    for (int i = 1; i < pop_size; i++) {
        if (a != population[i].points_fitted) {
            return 0;
        }
    }
    return 1;

}

/**
 * Escribe a configuración do plano atopada no ficheiro de saída "solutions.txt".
 * 
 * @param s Estructura solution do plano atopado.
 * 
 */
void write_solution(Solution s) {

    std::ofstream myfile;

    printf("Solución atopada: (%f,%f,%f,%f)\n", s.chromosome[0], s.chromosome[1],
            s.chromosome[2], s.chromosome[3]);
    printf("Puntos do plano: %i\n", s.points_fitted);
    printf("Escribindo resultados.\n");
    myfile.open("solutions.txt", std::ios::out | std::ios::app);
    myfile << "Solution : (" << s.chromosome[0] << "," << s.chromosome[1] << ","
            << s.chromosome[2] << "," << s.chromosome[3] << ")" << std::endl;
    myfile << "Points fitted: " << s.points_fitted << std::endl;
    myfile << "Points in region: " << s.points_in_region << std::endl;
    myfile.close();

}