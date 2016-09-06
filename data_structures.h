#ifndef DATA_STRUCTURES
#define DATA_STRUCTURES

//Estruturas xenéricas para puntos e solucións.

typedef struct Solution {
    float chromosome[4];
    int points_fitted;
    int points_in_region;
    float fitness;
} Solution;

typedef struct Point {
    float coordinates[3];  
} Point;

//Estrutura auxiliar para facer a ordenación.
struct SOLUTIONCmp {
    bool operator()(const Solution& s1, const Solution& s2) {
        return s2.fitness < s1.fitness;
    }
};

#endif
