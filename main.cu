
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <float.h>
#include <thrust/remove.h>
#include <thrust/execution_policy.h>
#include <fstream>
#include "data_structures.h"
#include "reduce.h"
#include "auxiliares.h"
#include "mutation.h"
#include "crossover.h"
#include "fitness.h"
#include "selection.h"
#include "generation.h"
#include "replacement.h"
#include "debug.h"

//SHUFFLE
#include <algorithm>    // std::shuffle
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock

#define THRESHOLD 0.012f 
#define PLANE_REGION 0.5f 
#define EPSILON FLT_EPSILON 
#define POPULATION_SIZE 200
#define CLOUD_SIZE 352600
#define MUTATION_D_INDEX 100 
#define MUTATION_RATE 0.25f
#define CROSS_D_INDEX 2.0f
#define CROSS_PROB 0.5f
#define GENE_SIZE 4.0f
#define MAX_EXECUTIONS 5000
#define SIGNIFICANT_PLANE 50
#define TOUR_SIZE 5
#define MAX_FAIL 3

/**
 * Función que executa o bucle principal do algoritmo. Ao contrario dos 
 * algoritmos xenéticos clásicos non atopa unha solución, senón que vai atopando
 * unha "parte" da solución . O criterio de parada do algoritmo é que todos os 
 * puntos da nube estén asignados a un plano (ou non queden máis planos) . O 
 * criterio de converxencia para a poboación é que todos os individuos da
 * poboación teñan o mesmo número de puntos axustados dentro deles.  
 * 
 * @param pop_size Tamaño da poboación.
 * @param cloud_size Tamaño da nube.
 * @param max_fail Número máximo de veces que se permite converxer ao algoritmo 
 * sen atopar un plano significativo.
 * @param population Punteiro a host_vector<Solution> que contén a poboación.
 * @param mating_pool Punteiro a host_vector<Solution> que contén a piscina de 
 * reprodución.
 * @param h_cloud Punteiro a host_vector<Point> que contén a nube en memoria 
 * principal do sistema.
 * @param d_cloud Punteiro a device_vector<Point> que contén a nube en memoria 
 * global de dispositivo.
 * @param plane_min_size Tamaño mínimo dun plano para consideralo 
 * arquitectónicamente significativo.
 * @param upper Punteiro ao array de límites superiores.
 * @param lower Punteiro ao array de límites inferiores.
 * @return 0 se o bucle se executa correctamente, o código de erro noutro caso.
 */
int bucle(size_t pop_size, size_t *cloud_size, int max_fail,
        thrust::host_vector<Solution> *population,
        thrust::host_vector<Solution> *mating_pool,
        thrust::host_vector<Point> *h_cloud,
        thrust::device_vector <Point> *d_cloud, size_t plane_min_size,
        float *upper, float *lower) {

    int fail = 0;
    int exec_count = 0;
    int result;
    float *mut_rate = (float*) malloc(sizeof (float));
    if (mut_rate == NULL) {
        printf("Non hai memoria suficiente para facer malloc\n");
        printf("Erro dentro do bucle\n");
        return (-1);
    }
    float *mut_d_index = (float*) malloc(sizeof (float));
    if (mut_d_index == NULL) {
        printf("Non hai memoria suficiente para facer malloc\n");
        printf("Erro dentro do bucle\n");
        return (-1);
    }
    Point *p_d_cloud = thrust::raw_pointer_cast(d_cloud->data());
    thrust::device_vector<int> fitted(*cloud_size);
    thrust::device_vector<int> region(*cloud_size);
    thrust::device_vector<float>fit(*cloud_size);

    //Criterio de parada do algoritmo
    while (*cloud_size > 0 && fail < max_fail) {
        *mut_rate = MUTATION_RATE;
        *mut_d_index = MUTATION_D_INDEX;
        exec_count = 0;
        thrust::fill(thrust::device, fitted.begin(), fitted.end(), 0);
        thrust::fill(thrust::device, region.begin(), region.end(), 0);
        thrust::fill(thrust::device, fit.begin(), fit.end(), 0.f);
        if ((result = generate_population(pop_size, *cloud_size, population,
                h_cloud)) != 0) //generation.cu
            return result;
        normalize(population, pop_size); //generation.cu
        evaluate_population_cuda(THRESHOLD, PLANE_REGION, EPSILON,
                p_d_cloud, population, *cloud_size, pop_size, &fitted, &region,
                &fit); //fitness.cu
        //Criterio de converxencia.
        while (!is_converged(*population, pop_size)) {//auxiliares.cu
            if ((result = tournament_selection(*population, mating_pool,
                    pop_size, TOUR_SIZE)) != 0) //selection.cu
                return -3;
            if ((result = recombine(mating_pool, pop_size, CROSS_PROB,
                    CROSS_D_INDEX, upper, lower, EPSILON)) != 0) //crossover.cu
                return -4;
            mutation(pop_size, mating_pool, mut_rate, mut_d_index, upper,
                    lower); //mutation.cu 
            evaluate_population_cuda(THRESHOLD, PLANE_REGION, EPSILON,
                    p_d_cloud, mating_pool, *cloud_size, pop_size, &fitted,
                    &region, &fit); //fitness.cu
            replacement(population, *mating_pool, pop_size); //replacement.cu
            exec_count++;
            recalculate_parameters(mut_d_index, mut_rate, GENE_SIZE, exec_count,
                    MAX_EXECUTIONS, MUTATION_D_INDEX); //mutation.cu
        }

        //Unha vez converxeu a poboación, retiramos da nube os puntos pertencentes
        //ao plano atopado.
        if (population->operator[](0).points_fitted > plane_min_size) {
            evaluate_population_cuda(THRESHOLD, PLANE_REGION, EPSILON, p_d_cloud,
                    population, *cloud_size, 1, &fitted, &region, 
                    &fit); //fitness.cu
            thrust::detail::normal_iterator<thrust::device_ptr<Point> > new_end
                    = thrust::remove_if(thrust::device, d_cloud->begin(),
                    d_cloud->end(), fitted.begin(), thrust::identity<int>());
            write_solution(population->operator[](0)); //auxiliares.cu
            *cloud_size = *cloud_size - population->operator[](0).points_fitted;
            fail = 0;
            thrust::copy(d_cloud->begin(), new_end, h_cloud->begin());
            exec_count = 0;
        } else {
            fail++;
        }
    }
    return 0;

}

/**
 * Función auxiliar que crea a nube de puntos sintética. Engade un delta 
 * aleatorio a cada punto para xerar ruído na nube.
 * 
 * @param cloud Punteiro a host_vector<Point> que contén a nube na memoria do 
 * sistema.
 * @param cloud_size Tamaño da nube
 * @param s  Estructura Solution cos parámetros do plano que vamos a empregar 
 * para xerar puntos
 * @param size Número de puntos a xerar.
 * @param offset Punto da nube onde empezamos a escribir.
 */
void generate_cloud(thrust::host_vector<Point> *cloud,Solution s, size_t size,
        size_t offset) {

    float z = 0;
    float delta = 0;
    for (int i = 0; i < 100; i++) {
        for (int j = 0; j < size / 100; j++) {
            z = ((-s.chromosome[3])+(-s.chromosome[0] * i)+
                    (-s.chromosome[1] * j)) / (s.chromosome[2]);
            delta = ((2 * ((float) rand() / (float) RAND_MAX)) - 1) / 200.0f;

            cloud->operator[](i * size / 100 + j + offset).coordinates[0] =
                    (float) i + delta;
            cloud->operator[](i * size / 100 + j + offset).coordinates[1] = 
                    (float) j + delta;
            cloud->operator[](i * size / 100 + j + offset).coordinates[2] = 
                    (float) z + delta;
        }
    }
    
}

int main() {
    srand(time(NULL));
    size_t * cloud_size = (size_t *) malloc(sizeof (size_t));
    *cloud_size = CLOUD_SIZE;
    thrust::host_vector<Point> h_cloud(*cloud_size);
    thrust::host_vector<Solution> population(POPULATION_SIZE);
    thrust::host_vector<Solution> mating_pool(POPULATION_SIZE);
    float *upper = (float *) malloc(4 * sizeof (float));
    float *lower = (float *) malloc(4 * sizeof (float));
    upper[0] = 1.0f;
    upper[1] = 1.0f;
    upper[2] = 1.0f;
    upper[3] = 1.0f;
    lower[0] = -1.0f;
    lower[1] = -1.0f;
    lower[2] = -1.0f;
    lower[3] = -1.0f;
    Solution s;

    s.chromosome[0] = 11;
    s.chromosome[1] = 16;
    s.chromosome[2] = 14;
    s.chromosome[3] = -15;
    size_t size = 20000;
    size_t offset = 0;
    generate_cloud(&h_cloud, s, size, offset);
    offset += size;

    s.chromosome[0] = 8;
    s.chromosome[1] = 2;
    s.chromosome[2] = 3;
    s.chromosome[3] = 1;
    size = 9500;
    generate_cloud(&h_cloud, s, size, offset);
    offset += size;

    s.chromosome[0] = -1;
    s.chromosome[1] = 4.5;
    s.chromosome[2] = 2.5;
    s.chromosome[3] = 12;
    size = 1200;
    generate_cloud(&h_cloud, s, size, offset);
    offset += size;

    s.chromosome[0] = 18;
    s.chromosome[1] = 22;
    s.chromosome[2] = 1;
    s.chromosome[3] = 9;
    size = 10500;
    generate_cloud(&h_cloud, s, size, offset);
    offset += size;

    s.chromosome[0] = -8;
    s.chromosome[1] = 0;
    s.chromosome[2] = 4;
    s.chromosome[3] = -3.25;
    size = 20000;
    generate_cloud(&h_cloud, s, size, offset);
    offset += size;

    s.chromosome[0] = -5.5;
    s.chromosome[1] = 12;
    s.chromosome[2] = -1;
    s.chromosome[3] = -1;
    size = 20000;
    generate_cloud(&h_cloud, s, size, offset);
    offset += size;

    s.chromosome[0] = 1;
    s.chromosome[1] = 1;
    s.chromosome[2] = 4;
    s.chromosome[3] = -5;
    size = 45000;
    generate_cloud(&h_cloud, s, size, offset);
    offset += size;

    s.chromosome[0] = 12;
    s.chromosome[1] = 3;
    s.chromosome[2] = 5.5;
    s.chromosome[3] = -1;
    size = 500;
    generate_cloud(&h_cloud,s, size, offset);
    offset += size;

    s.chromosome[0] = -1;
    s.chromosome[1] = 8;
    s.chromosome[2] = 2;
    s.chromosome[3] = -9;
    size = 2200;
    generate_cloud(&h_cloud, s, size, offset);
    offset += size;

    s.chromosome[0] = -3;
    s.chromosome[1] = -1;
    s.chromosome[2] = 8;
    s.chromosome[3] = -4;
    size = *cloud_size - offset;
    generate_cloud(&h_cloud, s, size, offset);
    offset += size;

    //Baraxamos aleatoriamente a poboación para enfrontar ao algoritmo ao peor 
    //caso posible (incrementa a diverxencia dentro dos warp no kernel de fitness)
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    shuffle(h_cloud.begin(), h_cloud.end(), std::default_random_engine(seed));
    thrust::device_vector<Point> d_cloud = h_cloud;
    int result=bucle(POPULATION_SIZE, cloud_size, MAX_FAIL, &population,
            &mating_pool, &h_cloud, &d_cloud, SIGNIFICANT_PLANE, upper, lower);
    if (result==1){
        cudaDeviceReset();
        exit(0);
    }
    cudaDeviceReset();
    exit (result);
}
