
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
#include <iostream>

//SHUFFLE
#include <algorithm>    // std::shuffle
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock

#define THRESHOLD 0.012f 
#define PLANE_REGION 0.5f 
#define EPSILON FLT_EPSILON 
#define POPULATION_SIZE 200
#define CLOUD_SIZE 100000
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
 * @param t Umbral do plano.
 * @param r Rexión do plano.
 * @param max_exec Número máximo de execucións permitida.  
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
 * @param tour_size Tamaño do torneo.
 * @param mutation_rate Probabilidade de mutación.
 * @param mutation_index Índice da función de distribución de probabilidade 
 * polinomial do operador de mutación.
 * @param cross_rate Probabilidade de cruce.
 * @param cross_index Índice da función de distribución de probabilidade do 
 * operador de cruce.
 * @param fast_convergence Booleano para indicar se se desexa activar a mellora
 * de converxencia anticipada.
 * @param min_growth Crecemento mínimo relativo aceptable no modo de 
 * converxencia anticipada. 
 * @param max_block_size Tamaño máximo dos bloques para os kernels de CUDA.
 * @param max_grid_size Tamaño máximo das mallas para os kernels de CUDA.
 * @return 0 se o algorimto se executa correctamente, código de erro noutro 
 * caso.
 */
int bucle(size_t pop_size, size_t *cloud_size, float t, float r,
        int max_exec, int max_fail, thrust::host_vector<Solution> *population,
        thrust::host_vector<Solution> *mating_pool,
        thrust::host_vector<Point> *h_cloud,
        thrust::device_vector <Point> *d_cloud, size_t plane_min_size,
        float *upper, float *lower, int tour_size, float mutation_rate,
        float mutation_index, float cross_rate, float cross_index,
        bool fast_convergence, float min_growth, int max_block_size,
        int max_grid_size) {

    int planes_found = 1;
    float av_fitness;
    float time_acumulated = 0;
    std::ofstream myfile;
    myfile.open("performance.txt", std::ios::out | std::ios::trunc);
    cudaEvent_t start, stop, eval_start, eval_stop;
    float miliseconds;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventCreate(&eval_start);
    cudaEventCreate(&eval_stop);
    cudaDeviceProp prop;
    int device;
    gpuErrchk(cudaGetDevice(&device));
    gpuErrchk(cudaGetDeviceProperties(&prop, device));

    //cambio de criterio de convergencia
    float * average = (float*) malloc(sizeof (float));
    *average = 1;
    int * not_improved = (int *) malloc(sizeof (int));

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
    myfile << "Inicio do algoritmo." << std::endl;
    myfile << "Tamaño da nube: " << *cloud_size << std::endl;
    myfile << "Tamaño da poboación: " << pop_size << std::endl;
    myfile << "Tamaño máximo de bloque (CUDA): " << max_block_size << std::endl;
    myfile << "Tamaño máximo de malla (CUDA): " << max_grid_size << std::endl;
    //Criterio de parada do algoritmo
    while (*cloud_size > 0 && fail < max_fail) {
        cudaEventRecord(start);
        *mut_rate = mutation_rate;
        *mut_d_index = mutation_index;
        exec_count = 0;
        thrust::fill(thrust::device, fitted.begin(), fitted.end(), 0);
        thrust::fill(thrust::device, region.begin(), region.end(), 0);
        thrust::fill(thrust::device, fit.begin(), fit.end(), 0.f);
        if ((result = generate_population(pop_size, *cloud_size, population,
                h_cloud)) != 0) //generation.cu
            return result;
        myfile << "#################################" << std::endl;
        myfile << "Plano número " << planes_found << std::endl;
        normalize(population, pop_size); //generation.cu
        evaluate_population_cuda(t, r, EPSILON,
                p_d_cloud, population, *cloud_size, pop_size, &fitted, &region,
                &fit, max_block_size, max_grid_size, device, prop); //fitness.cu
        //Criterio de converxencia.
        while (!is_converged(*population, pop_size, average, not_improved,
                fast_convergence, min_growth)) {//auxiliares.cu
            if ((result = tournament_selection(*population, mating_pool,
                    pop_size, tour_size)) != 0) //selection.cu
                return -3;
            if ((result = recombine(mating_pool, pop_size, cross_rate,
                    cross_index, upper, lower, EPSILON)) != 0) //crossover.cu
                return -4;
            mutation(pop_size, mating_pool, mut_rate, mut_d_index, upper,
                    lower); //mutation.cu 
            cudaEventRecord(eval_start);
            evaluate_population_cuda(t, r, EPSILON, p_d_cloud, mating_pool,
                    *cloud_size, pop_size, &fitted, &region, &fit,
                    max_block_size, max_grid_size, device, prop); //fitness.cu
            cudaEventRecord(eval_stop);
            cudaEventSynchronize(eval_stop);
            cudaEventElapsedTime(&miliseconds, eval_start, eval_stop);
            time_acumulated += miliseconds;
            myfile << "Tempo para avaliar a poboación, ciclo " << exec_count
                    << ": " << miliseconds << "ms" << std::endl;
            replacement(population, *mating_pool, pop_size); //replacement.cu
            av_fitness = average_fitness(*population, pop_size);
            myfile << "Fitness medio da poboación: " << av_fitness << std::endl;
            exec_count++;
            recalculate_parameters(mut_d_index, mut_rate, GENE_SIZE, exec_count,
                    max_exec, mutation_index); //mutation.cu
        }

        //Unha vez converxeu a poboación, retiramos da nube os puntos
        //pertencentes ao plano atopado.
        if (population->operator[](0).points_fitted > plane_min_size) {
            evaluate_population_cuda(t, r, EPSILON, p_d_cloud,
                    population, *cloud_size, 1, &fitted, &region,
                    &fit, max_block_size, max_grid_size, device,
                    prop); //fitness.cu
            thrust::detail::normal_iterator<thrust::device_ptr<Point> > new_end
                    = thrust::remove_if(thrust::device, d_cloud->begin(),
                    d_cloud->end(), fitted.begin(), thrust::identity<int>());
            cudaEventRecord(stop);
            cudaEventSynchronize(stop);
            cudaEventElapsedTime(&miliseconds, start, stop);
            myfile << "Tempo de execución total para atopar o plano "
                    << planes_found << ": " << miliseconds << " ms" <<
                    std::endl;
            myfile << "Tempo de execución de avaliación de fitnes para este"
                    " plano: " << time_acumulated << " ms ." << std::endl;
            myfile << "Tempo medio de avaliación por ciclo: " <<
                    time_acumulated / exec_count << " ms." << std::endl;
            myfile << "Porcentaxe de tempo respecto do tempo total: " <<
                    time_acumulated / miliseconds * 100 << " %" << std::endl;
            myfile << "#################################" << std::endl;
            time_acumulated = 0;
            write_solution(population->operator[](0)); //auxiliares.cu
            *cloud_size = *cloud_size - population->operator[](0).points_fitted;
            fail = 0;
            thrust::copy(d_cloud->begin(), new_end, h_cloud->begin());
            exec_count = 0;
            planes_found++;

        } else {
            time_acumulated=0;
            exec_count=0;
            fail++;
        }
    }
    myfile.close();
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
void generate_cloud(thrust::host_vector<Point> *cloud, Solution s, size_t size,
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

int main(int argc, char **argv) {
    srand(time(NULL));
    size_t * cloud_size = (size_t *) malloc(sizeof (size_t));
    size_t pop_size;
    size_t offset = 0;
    size_t size;
    float mut_rate;
    float mut_d_index;
    float cross_d_index;
    float cross_rate;
    float r;
    float t;
    int max_exec;
    int significant_plane;
    int tour_size;
    int max_fail;
    bool fast_convergence;
    float min_growth = 0.f;
    float *upper = (float *) malloc(4 * sizeof (float));
    float *lower = (float *) malloc(4 * sizeof (float));
    int max_block_size;
    int max_grid_size;
    upper[0] = 1.0f;
    upper[1] = 1.0f;
    upper[2] = 1.0f;
    upper[3] = 1.0f;
    lower[0] = -1.0f;
    lower[1] = -1.0f;
    lower[2] = -1.0f;
    lower[3] = -1.0f;
    Solution s;


    switch (argc) {
        case 1: std::cout << "Modo normal" << std::endl;
            *cloud_size = CLOUD_SIZE;
            pop_size = POPULATION_SIZE;
            mut_rate = MUTATION_RATE;
            mut_d_index = MUTATION_D_INDEX;
            cross_d_index = CROSS_D_INDEX;
            cross_rate = CROSS_PROB;
            r = PLANE_REGION;
            t = THRESHOLD;
            max_exec = MAX_EXECUTIONS;
            significant_plane = SIGNIFICANT_PLANE;
            tour_size = TOUR_SIZE;
            max_fail = MAX_FAIL;
            fast_convergence = false;
            max_block_size = 256;
            max_grid_size = 128;
            break;
        case 7: std::cout << "Modo rendemento" << std::endl;
            if ((pop_size = std::stoi(argv[1])) % 2) {
                std::cout << "A poboación non pode ser impar" << std::endl;
                exit(-1);
            }
            if ((*cloud_size = std::stoi(argv[2])) % 10000) {
                std::cout << "O tamaño da nube de puntos debe ser múltiplo"
                        " de 10000" << std::endl;
                exit(-1);
            }
            if (std::stoi(argv[3])) {
                fast_convergence = true;
            } else
                fast_convergence = false;
            min_growth = std::stof(argv[4]);
            max_block_size = std::stoi(argv[5]);
            if (!isPow2(max_block_size) || max_block_size > 1024) {
                std::cout << "Tamaño de bloque non válido!!" << std::endl;
                exit(-1);
            }
            max_grid_size = std::stoi(argv[6]);
            if (!isPow2(max_grid_size) || max_grid_size > 65535) {
                std::cout << "Tamaño de malla non válido!!" << std::endl;
                exit(-1);
            }
            mut_rate = MUTATION_RATE;
            mut_d_index = MUTATION_D_INDEX;
            cross_d_index = CROSS_D_INDEX;
            cross_rate = CROSS_PROB;
            r = PLANE_REGION;
            t = THRESHOLD;
            max_exec = MAX_EXECUTIONS;
            significant_plane = SIGNIFICANT_PLANE;
            tour_size = TOUR_SIZE;
            max_fail = MAX_FAIL;
            break;
        case 17: std::cout << "Modo avanzado" << std::endl;
            if ((pop_size = std::stoi(argv[1])) % 2) {
                std::cout << "A poboación non pode ser impar" << std::endl;
                exit(-1);
            }
            if ((*cloud_size = std::stoi(argv[2])) % 10000) {
                std::cout << "O tamaño da nube de puntos debe ser múltiplo "
                        "de 10000" << std::endl;
                exit(-1);
            }
            r = std::stof(argv[3]);
            t = std::stof(argv[4]);
            if (r < t) {
                std::cout << "O parámetro r non pode ser menor que o parámetro "
                        "t" << std::endl;
                exit(-1);
            }
            mut_d_index = std::stof(argv[5]);
            if ((mut_rate = std::stof(argv[6])) > 1.0f) {
                std::cout << "A probabilidade de mutación non pode ser maior "
                        "ca 1!!" << std::endl;
                exit(-1);
            }
            cross_d_index = std::stof(argv[7]);

            if ((cross_rate = std::stof(argv[8])) > 1.0f) {
                std::cout << "A probabilidade de cruzamento non pode ser maior"
                        " ca 1!!" << std::endl;
                exit(-1);
            }
            max_exec = std::stoi(argv[9]);
            significant_plane = std::stoi(argv[10]);
            tour_size = std::stoi(argv[11]);
            max_fail = std::stoi(argv[12]);
            if (std::stoi(argv[13])) {
                fast_convergence = true;
            } else
                fast_convergence = false;
            min_growth = std::stof(argv[14]);
            max_block_size = std::stoi(argv[15]);
            if (!isPow2(max_block_size) || max_block_size > 1024) {
                std::cout << "Tamaño de bloque non válido!!" << std::endl;
                exit(-1);
            }
            max_grid_size = std::stoi(argv[16]);
            if (!isPow2(max_grid_size) || max_grid_size > 65535) {
                std::cout << "Tamaño de malla non válido!!" << std::endl;
                exit(-1);
            }
            break;
        default: std::cout << "Número de parámetros incorrecto" << std::endl;
            std::cout << "Uso: " << std::endl;
            std::cout << "Modo normal: Execución cos parámetros preestablecidos"
                    << std::endl;
            std::cout << "./paralelo" << std::endl;
            std::cout << "##################################################"
                    "####"
                    << std::endl;
            std::cout << "Modo rendemento: Estuda o rendemento do algoritmo "
                    "cambiando os parámetros que afectan ao tempo de execución"
                    << std::endl;
            std::cout << "./paralelo pop_size cloud_size fast_convergence "
                    "min_growth max_block_size max_grid_size" << std::endl;
            std::cout << "\tpop_size: tamaño da poboación (número par)"
                    << std::endl;
            std::cout << "\tcloud_size: tamaño da nube (múltiplo de 10000)"
                    << std::endl;
            std::cout << "\tfast_convergence 1 para activar a converxencia "
                    "adiantada, 0 para desactivala" << std::endl;
            std::cout << "\tmin_growth crecemento mínimo aceptable . "
                    "Recomendado"
                    "[0.00001,0.001] . (Con fast_convergence=0, introducir "
                    "calquera valor." << std::endl;
            std::cout << "\tmax_block_size: Tamaño de bloque dos kernels "
                    "(múltiplo de 2, menor de 1024)" << std::endl;
            std::cout << "\tmax_grid_size: Tamaño máximo da malla (múltiplo de"
                    " 2, menor de 65535" << std::endl;
            std::cout << "############################################"
                    "##########"
                    << std::endl;
            std::cout << "Modo avanzado:" << std::endl;
            std::cout << "./paralelo pop_size cloud_size r t mut_d_index mut_rate"
                    "cross_d_index cross_rate max_exec significant_plane "
                    "tour_size max_fail fast_convergence min_growth "
                    "max_block_size max_grid_size" << std::endl;
            exit(-1);
            break;
    }

    thrust::host_vector<Point> h_cloud(*cloud_size);
    thrust::host_vector<Solution> population(pop_size);
    thrust::host_vector<Solution> mating_pool(pop_size);

    s.chromosome[0] = 11;
    s.chromosome[1] = 16;
    s.chromosome[2] = 14;
    s.chromosome[3] = -15;
    size = (*cloud_size / 100)*20;

    generate_cloud(&h_cloud, s, size, offset);
    offset += size;

    s.chromosome[0] = 8;
    s.chromosome[1] = 2;
    s.chromosome[2] = 3;
    s.chromosome[3] = 1;
    size = (*cloud_size / 100)*15;
    generate_cloud(&h_cloud, s, size, offset);
    offset += size;

    s.chromosome[0] = -1;
    s.chromosome[1] = 4.5;
    s.chromosome[2] = 2.5;
    s.chromosome[3] = 12;
    size = (*cloud_size / 100)*12;
    generate_cloud(&h_cloud, s, size, offset);
    offset += size;

    s.chromosome[0] = 18;
    s.chromosome[1] = 22;
    s.chromosome[2] = 1;
    s.chromosome[3] = 9;
    size = (*cloud_size / 100)*10;
    generate_cloud(&h_cloud, s, size, offset);
    offset += size;

    s.chromosome[0] = -8;
    s.chromosome[1] = 0;
    s.chromosome[2] = 4;
    s.chromosome[3] = -3.25;
    size = (*cloud_size / 100)*10;
    generate_cloud(&h_cloud, s, size, offset);
    offset += size;

    s.chromosome[0] = -5.5;
    s.chromosome[1] = 12;
    s.chromosome[2] = -1;
    s.chromosome[3] = -1;
    size = (*cloud_size / 100)*9;
    generate_cloud(&h_cloud, s, size, offset);
    offset += size;

    s.chromosome[0] = 1;
    s.chromosome[1] = 1;
    s.chromosome[2] = 4;
    s.chromosome[3] = -5;
    size = (*cloud_size / 100)*8;
    generate_cloud(&h_cloud, s, size, offset);
    offset += size;

    s.chromosome[0] = 12;
    s.chromosome[1] = 3;
    s.chromosome[2] = 5.5;
    s.chromosome[3] = -1;
    size = (*cloud_size / 100)*7;
    generate_cloud(&h_cloud, s, size, offset);
    offset += size;

    s.chromosome[0] = -1;
    s.chromosome[1] = 8;
    s.chromosome[2] = 2;
    s.chromosome[3] = -9;
    size = (*cloud_size / 100)*6;
    generate_cloud(&h_cloud, s, size, offset);
    offset += size;

    s.chromosome[0] = -3;
    s.chromosome[1] = -1;
    s.chromosome[2] = 8;
    s.chromosome[3] = -4;
    size = (*cloud_size / 100)*3;
    generate_cloud(&h_cloud, s, size, offset);
    offset += size;

    //Baraxamos aleatoriamente a poboación para enfrontar ao algoritmo ao peor 
    //caso posible (incrementa a diverxencia dentro dos warp no kernel 
    //de fitness)
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    shuffle(h_cloud.begin(), h_cloud.end(), std::default_random_engine(seed));
    thrust::device_vector<Point> d_cloud = h_cloud;


    int result = bucle(pop_size, cloud_size, t, r, max_exec, max_fail,
            &population, &mating_pool, &h_cloud, &d_cloud, significant_plane,
            upper, lower, tour_size, mut_rate, mut_d_index, cross_rate,
            cross_d_index, fast_convergence, min_growth, max_block_size,
            max_grid_size);

    if (result == 1) {
        cudaDeviceReset();
        exit(0);
    }
    cudaDeviceReset();
    exit(result);
}
