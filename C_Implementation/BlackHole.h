#ifndef __BH__ 
#define __BH__ 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct line Linha_de_mundo; 

/*
x0 = [t,r,theta,phi,pr,ptheta]
m = Massa da partícula
q = carga da partícula
E = energia da partícula
lz = Momento angular no eixo z da partícula
Qcarter = constante de Carter da partícula
M = massa do buraco negro
Q = carga do buraco negro
a = Momento angular por massa do buraco negro

*/
Linha_de_mundo* create_Wline(size_t N, double dt, double* x0, double m, double q, double E, double lz, double Qcarter, double M, double Q, double a);

Linha_de_mundo* reuse_Wline(Linha_de_mundo* ln,size_t N, double dt, double* x0, double m, double q, double E, double lz, double Qcarter, double M, double Q, double a);

// Evolui Linha_de_mundo 
double* evolve_RK(Linha_de_mundo* ln);

// Evolui Linha_de_mundo ao contrário
double* evolve_RK_inverse(Linha_de_mundo* ln);
// Libera Linha_de_mundo da memória
size_t freeWL(Linha_de_mundo* ln);

// Retorna posição final da partícula 
double* get_x(Linha_de_mundo* ln);

// Retorna posição inicial da partícula
double* get_x0(Linha_de_mundo* ln);

#endif