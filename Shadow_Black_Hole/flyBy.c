#include "BlackHole.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef PI
#define PI 3.14159265358979323846
#endif


int main(void){
    const double M = 1.0; 
    const double a = sqrt(0.84); 
    const double Q = 0.2;

    double E = sqrt(2.0); 
    const double Lz = 4.2;
    const double m = 1.0;
    const double q = 0.0;
    double x0[6] = {0.0,13.0, PI/2, 0.0,-1,1};
    size_t N = 10000;
    double interval = (E-0.9685)/N;
    for (size_t i = 0; i<(N+1); i++){
        Linha_de_mundo* ln = create_Wline(50000,0.1,x0,m,q,E,Lz,0,M,Q,a); 
        double* xf = evolve_RK(ln);
        printf("Posição Final %zu ",i);
        freeWL(ln);
        E = E - interval;
    }
}