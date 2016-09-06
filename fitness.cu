#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include "data_structures.h"
#include <float.h>
#include "auxiliares.h"
#include "reduce.h"

/**
 * Función auxiliar para código de dispositivo que calcula o identificador
 * global do thread en unha dimensión.
 * 
 * @return Un enteiro co identificador global do thread
 */
__device__ int getGlobalIdx_1D_1D() {
    return blockIdx.x * blockDim.x + threadIdx.x;
}

/**
 * Kernel cuda para o cálculo do fitness de un individuo da poboación. Obtén a 
 * distancia euclídea de un punto respecto do plano hipotetizado na solución. 
 * Calcula o fitness parcial de cada punto e establece se dito punto pertence á 
 * rexión ou ao umbral do plano.
 * 
 * @param t Umbral do plano (threshold).
 * @param r Rexión do plano (region).
 * @param eps Épsilon de máquina.
 * @param cloud Punteiro á nube en memoria global do dispositivo.
 * @param s Estrutura solution do plano a evaluar.
 * @param fitted Array para conteo de puntos pertencentes ao plano.
 * @param region Array para conteo de puntos na rexión do plano.
 * @param fit Array para almacenar os fitness parciais de cada punto.
 * @param b Base do logaritmo.
 * @param max_log Valor máximo da función logarítmica.
 */
__global__ void
fitness(float t, float r, float eps, Point *cloud, Solution *s, int * fitted,
        int * region, float * fit, float b, float max_log) {

    //Arrays en memoria compartida
    extern __shared__ float array[];
    float *d = (float*) array;
    float *k = (float*) &d[blockDim.x];
    int i = getGlobalIdx_1D_1D();
    int tid = threadIdx.x;
    //Cálculo da distancia euclídea.Facemos as escrituras en memoria compartida.
    k[tid] = sqrtf(powf(s->chromosome[0], 2.0f) + powf(s->chromosome[1], 2.0f)
            + powf(s->chromosome[2], 2.0f));
    d[tid] = fabsf(cloud[i].coordinates[0] * s->chromosome[0] + 
            cloud[i].coordinates[1]* s->chromosome[1] + cloud[i].coordinates[2]
            * s->chromosome[2] + s->chromosome[3]);
    if (k[tid] <= 5 * FLT_EPSILON)
        d[tid] = d[tid] / eps;
    else
        d[tid] = d[tid] / k[tid];
    //Facemos as asignacións empregando predicación de instruccións (guía de 
    //programación de cuda, punto 5.4.2). Esta opción é un 12.5 % máis rápida ca
    //empregar sentencias if-else e un 25% máis rápida que a versión sen
    //branching
    //(region[i]=1+r-d[tid], fit[i]=1+t-d[tid]);
    if (d[tid] > 5 * FLT_EPSILON) {
        region[i] = (d[tid] <= r ? 1 : 0);
        fit[i] = d[tid] <= r ? ((log10f(d[tid]) / log10f(b))-(log10f(t) / 
                log10f(b)) - max_log) : 0;
        fitted[i] = (d[tid] <= t ? 1 : 0);
    } else {
        region[i] = 1;
        fitted[i] = 0;
        fit[i] = 0.f;
    }

}

/**
 * Función auxiliar para calcular o tamaño de bloque para as chamadas ao kernel 
 * de evaluación de fitness.
 * @param original Tamaño actual do bloque.
 * @param max Tamaño máximo de bloque.
 * @return O novo tamaño de bloque.
 */
int getBlockSize(int original, int max) {

    while (original > max) {
        original /= 2;
    }
    return original;

}

/**
 * Función auxiliar que chama secuncialmente ao kernel de fitness. 
 * Posteriormente, chama ao kernel de reducción e actualiza os valores de 
 * fitness, puntos na rexión e puntos encaixados no plano.
 * 
 * @param t Umbral do plano (threshold).
 * @param r Rexión do plano (region).
 * @param eps Epsilon de máquina.
 * @param p_d_cloud Punteiro á nube en memoria global
 * @param population Punteiro a host_vector<Solution> que contén a poboación.
 * @param cloud_size Tamaño da nube.
 * @param pop_size Tamaño da poboación.
 * @param fitted Array para conteo de puntos pertencentes ao plano.
 * @param region Array para conteo de puntos na rexión do plano.
 * @param fit Array para almacenar os fitness parciais de cada punto.
 * @param maxblock Tamaño máximo do bloque para os kernels.
 * @param maxgrid Tamaño máxima da malla de bloques.  
 * @return 0 se a execución é correcta, ou código de erro noutro caso. 
 */
void evaluate_population_cuda(float t, float r, float eps, Point *p_d_cloud, 
        thrust::host_vector<Solution> *population, size_t cloud_size,
        size_t pop_size, thrust::device_vector <int> *fitted, 
        thrust::device_vector <int> *region, thrust::device_vector<float>*fit,
        int maxblock, int maxgrid) {

    int gridSize = 0;
    int blockSize = 0;
    int size = cloud_size;
    int s = cloud_size;
    size_t offset = 0;
    size_t shmemSize = 0;
    float b = powf((r / t), (1.0f / r));
    float max_log = (log10f(r) / log10f(b))-(log10f(t) / log10f(b));
    int *p_fitted = thrust::raw_pointer_cast(&fitted->operator[](0));
    int *p_region = thrust::raw_pointer_cast(&region->operator[](0));
    float *p_fit = thrust::raw_pointer_cast(&fit->operator[](0));
    thrust::device_vector<Solution> d_pop = *population;
    Solution *p_d_pop = thrust::raw_pointer_cast(&d_pop.operator[](0));
    getNumBlocksAndThreads(size, maxgrid, maxblock, gridSize, 
            blockSize);
    thrust::host_vector<int> h_salida(2 * pop_size * sizeof (int));
    thrust::host_vector<float>h_fitness(pop_size * sizeof (float));
    thrust::device_vector<int> d_fitted_salida(gridSize * (sizeof (int)));
    int *p_d_fitted_salida = thrust::raw_pointer_cast(&d_fitted_salida[0]);
    thrust::device_vector<int> d_region_salida(gridSize * (sizeof (int)));
    int *p_d_region_salida = thrust::raw_pointer_cast(&d_region_salida[0]);
    thrust::device_vector<float> d_fit_salida(gridSize * (sizeof (int)));
    float *p_d_fit_salida = thrust::raw_pointer_cast(&d_fit_salida[0]);
    thrust::device_vector <int> salida_datos_enteros(2 * pop_size 
    * sizeof (int));
    int *p_salida_datos_enteros = thrust::raw_pointer_cast
    (&salida_datos_enteros[0]);
    thrust::device_vector <float>salida_datos_float(pop_size * sizeof (int));
    float *p_salida_datos_float = thrust::raw_pointer_cast
    (&salida_datos_float[0]);

    for (int j = 0; j < pop_size; j++) {
        s = cloud_size;
        offset = 0;
        blockSize = maxblock;
        while (s >= 1) {
            shmemSize = 2 * blockSize * sizeof (float);
            gridSize = s / blockSize;
            fitness << <gridSize, blockSize, shmemSize >> >(t, r, eps, p_d_cloud
                    + offset, p_d_pop + j, p_fitted + offset, p_region + offset,
                    p_fit + offset, b, max_log);
            gpuErrchk(cudaPeekAtLastError());
            offset += blockSize*gridSize;
            s = s % blockSize;
            blockSize = getBlockSize(blockSize, s);
        }
        getNumBlocksAndThreads(size, maxgrid, maxblock, gridSize,
                blockSize);
        performReduction(size, blockSize, gridSize, maxblock, maxgrid,
                p_salida_datos_enteros, p_fitted, p_d_fitted_salida, 2 * j);
        performReduction(size, blockSize, gridSize, maxblock, maxgrid,
                p_salida_datos_enteros, p_region, p_d_region_salida, 2 * j + 1);
        performReduction(size, blockSize, gridSize, maxblock, maxgrid,
                p_salida_datos_float, p_fit, p_d_fit_salida, j);
    }
    thrust::copy(salida_datos_enteros.begin(), salida_datos_enteros.end(), 
            h_salida.begin());
    thrust::copy(salida_datos_float.begin(), salida_datos_float.end(),
            h_fitness.begin());
    for (int j = 0; j < pop_size; j++) {
        population->operator[](j).points_fitted = h_salida[2 * j];
        population->operator[](j).points_in_region = h_salida[2 * j + 1];
        if (population->operator[](j).points_in_region == 0)
            population->operator[](j).fitness = -INFINITY;
        else
            population->operator[](j).fitness = (-(h_fitness[j])) /
                    population->operator[](j).points_in_region *
                    population->operator[](j).points_fitted;
    }

}
