#ifndef FITNESS_H
#define FITNESS_H

float fitness(float t, float r, float eps, std::vector<Point> cloud,
        float b, Solution * p, size_t cloud_size);

void evaluate_population(float t, float r, float eps,
        std::vector<Point> cloud, float b,
        std::vector<Solution> *population, size_t cloud_size,
        size_t pop_size) ;


#endif