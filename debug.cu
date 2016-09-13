#include <stdio.h>
#include "data_structures.h"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <fstream>
#include <thrust/count.h>
#include <thrust/execution_policy.h>



/*
 Auxiliar
 * Muestra un cromosoma
 */
void show_chromosome(Solution * s) {

    for (int i = 0; i < 4; i++) {
        printf("%f ", s->chromosome[i]);
    }
    printf("\n");

}

/*
 Auxiliar
 *Muestra una población
 */

void show_population(thrust::host_vector<Solution>genome, size_t pop_size) {

    for (int i = 0; i < pop_size; i++) {
        printf("******************\n");
        printf("Chromosoma [%i]: ", i);
        printf("(%f,%f,%f,%f)\n", genome[i].chromosome[0], 
                genome[i].chromosome[1], genome[i].chromosome[2], 
                genome[i].chromosome[3]);
        printf("Points fitted: %i\n", genome[i].points_fitted);
        printf("Points in region:%i\n",genome[i].points_in_region);
        printf("With fitness: %f\n", genome[i].fitness);
//        printf("ADN: %i%i%i\n",genome[i].adn[0],genome[i].adn[1],
        genome[i].adn[2]);
    }
    printf("******************\n");
}

/**
 * Auxiliar
 * Muestra la nube de puntos
 * @param cloud
 * @param cloud_size
 */
void show_cloud(thrust::host_vector<Point> cloud, size_t cloud_size) {

    for (int i = 0; i < cloud_size; i++) {
        printf("Point[%i]: (%f,%f,%f)\n", i, cloud[i].coordinates[0], 
                cloud[i].coordinates[1], cloud[i].coordinates[2]);
    }
}

void write_eliminated(thrust::host_vector<Point> cloud, size_t cloud_size, 
        const char *str, thrust::device_vector<int> vector
/*,std::ios::openmode mode*/) {
    thrust::host_vector<int> eliminados = vector;

    std::ofstream myfile;
    //mode=std::ios::out|std::ios::ate;
    myfile.open(str, std::ios::out | std::ios::app);
    for (int i = 0; i < cloud_size; i++) {
        if (eliminados[i] == 1) {
            myfile << "Point:[" << i << "]:(" << cloud[i].coordinates[0] 
                    << "," << cloud[i].coordinates[1] << "," 
                    << cloud[i].coordinates[2] << ")" << std::endl;
            //myfile << "Origin: " << cloud[i].origin << std::endl;

        }
    }
    myfile.close();


}
void write_eliminated_serie(thrust::host_vector<Point> cloud, size_t eliminados,
        const char *str /*,std::ios::openmode mode*/) {

    size_t size=cloud.size();
    std::ofstream myfile;
    //mode=std::ios::out|std::ios::ate;
    myfile.open(str, std::ios::out | std::ios::app);
    for (int i = size - eliminados; i < size; i++) {
        
            myfile << "Point:[" << i << "]:(" << cloud[i].coordinates[0] << "," 
                    << cloud[i].coordinates[1] << "," << 
                    cloud[i].coordinates[2] << ")" << std::endl;
            //myfile << "Origin: " << cloud[i].origin << std::endl;

    }
    myfile.close();


}


void write_cloud(thrust::host_vector<Point> cloud, size_t cloud_size,
        const char *str/*,std::ios::openmode mode*/) {

    std::ofstream myfile;
    //mode=std::ios::out|std::ios::ate;
    myfile.open(str, std::ios::out | std::ios::trunc);
    for (int i = 0; i < cloud_size; i++) {
        myfile << "Point:[" << i << "]:(" << cloud[i].coordinates[0] << ","
                << cloud[i].coordinates[1] << ","
                << cloud[i].coordinates[2] << ")" << std::endl;
   //     myfile << "Origin: " << cloud[i].origin << std::endl;

    }
    myfile.close();


}

void write_population(thrust::host_vector<Solution> population, size_t pop_size,
        const char *str/*,std::ios::openmode mode*/) {

    std::ofstream myfile;
    //mode=std::ios::out|std::ios::ate;
    myfile.open(str, std::ios::out | std::ios::trunc);
    for (int i = 0; i < pop_size; i++) {
        myfile << "Solution : (" << population[i].chromosome[0] << "," 
                << population[i].chromosome[1] << "," 
                << population[i].chromosome[2] << "," 
                << population[i].chromosome[3] << ")" << std::endl;
        myfile << "Fitness: "<<population[i].fitness<<std::endl;
        myfile << "Points fitted: " <<
                population[i].points_fitted << std::endl;
        myfile << "Points in region: "
                <<population[i].points_in_region<<std::endl;
     //   myfile << "Origin: "<<population[i].adn[0]<<population[i].adn[1]
        <<population[i].adn[2]<<std::endl;


    }
    myfile.close();
}

int contar(thrust::device_vector<int> vector) {


    return (thrust::count(thrust::device, vector.begin(), vector.end(), 1));

}

void show_vector(thrust::device_vector <int> vector, size_t size) {
    thrust::host_vector<int> copia = vector;

    for (int i = 0; i < size; i++) {
        printf("%i\n", copia[i]);

    }
}

void write_vector (thrust::device_vector<int>vector,size_t size,
        const char *str){

    thrust::host_vector<int> v=vector;
 std::ofstream myfile;
    //mode=std::ios::out|std::ios::ate;
    myfile.open(str, std::ios::out | std::ios::ate);
    for (int i = 0; i < size; i++) {
        myfile<< v[i]<<std::endl;
    }
    myfile.close();

}