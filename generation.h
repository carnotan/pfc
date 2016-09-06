#ifndef GENERATION_H
#define GENERATION_H

void get_solution(Point p1, Point p2, Point p3, Solution *s);
void normalize(std::vector<Solution> *population, size_t pop_size);
int generate_population(size_t pop_size, size_t cloud_size, 
        std::vector<Solution> *population, 
        std::vector<Point> *cloud);

#endif