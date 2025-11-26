#include "DBlackHole.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef PI
#define PI 3.14159265358979323846
#endif


int main(void){
    const double M = 1.0; 
    //const double a = sqrt(0.84); 
    //const double a = 0.0;
    const double a = 0.99;
    const double Q = 0.0;

    // Partícula tipo-luz
    const double m = 0.0;
    const double q = 0.0;

    // Posição do observador
    const double r0 = 200;
    //const double theta0 = 0.1;
    //const double theta0 = 0;
    const double theta0 =  0.01 ;
    printf("theta %lf \n",theta0);

    // Termos auxiliares na conta
    const double delta = r0*r0 + a*a + Q*Q - 2*M*r0;
    const double sigma = r0*r0 + a*a*cos(theta0)*cos(theta0);
    // Métrica
    const double g_tt = - (delta - a*a*sin(theta0)*sin(theta0))/sigma;
    const double g_tphi = -a*sin(theta0)*sin(theta0)*(r0*r0 + a*a - delta)/sigma;
    const double g_phiphi = ((r0*r0+a*a)*(r0*r0+a*a) - delta*a*a*sin(theta0)*sin(theta0))*sin(theta0)*sin(theta0)/sigma;
    const double g_rr = sigma/delta;
    const double g_thetatheta = sigma;

    // Parãmetros da simulação
    const size_t Nx = 300;
    const size_t Ny = 300;
    const double x_interval = 20.0; // -10 -> 10
    const double y_interval = 20.0; // 10 -> 10

    // Termos repetitivos nos LOOPS
    const double dx = (2.0*x_interval)/Nx;
    const double dy = (2.0*y_interval)/Ny;
    const double termx = 1.0/(sqrt(g_phiphi)*r0);
    const double termy = 1.0/(sqrt(g_thetatheta)*r0);
    // raio do buraco negro
    const double r_bh = M + sqrt(M*M -a*a -Q*Q);

    printf("Raio do Buraco Negro %lf \n",r_bh);



    // Número de trajetórias
    size_t N = Nx*Ny;
    double x0[6] = {0.0,r0,theta0,0.0,0.0,0.0};
    double* shadow = (double*) malloc(N*sizeof(double));
    Linha_de_mundo* ln = create_Wline(100000,0.001,x0,0.0,0.0,0.0,0.0,NAN,0.0,0.0,0.0); 
    // Recursão em y
    for (size_t j = 0; j<Ny; j++){
        // Recursão em x
        for(size_t i = 0; i<Nx; i++){
            //printf("Iniciando trajetória %zu %zu \n",j,i);
            // g_t,phi = -2*a(r**2 + a**2)
            double x = i*dx - x_interval;
            double y = j*dy - y_interval;
            double Lz = -x*termx*g_phiphi + g_tphi; 
            double E = -(-g_tphi*termx*x + g_tt);
            double ptheta = y*sqrt(g_thetatheta)/r0;
            double pr = sqrt((-g_tt - g_phiphi*termx*termx*x*x - 2*g_tphi*termx*x - g_thetatheta*termy*termy*y*y)*g_rr);

            //printf("pr %lf\n",pr);
            x0[4] = pr;
            x0[5] = ptheta;
            ln = reuse_Wline(ln,50000,0.01,x0,m,q,E,Lz,NAN,M,Q,a); 
            double* xf = evolve_RK_inverse(ln);
            if((xf[1]<r_bh*1.0001)||isnan(xf[1])){
                shadow[Nx*j + i] = 0;
                //printf("Trajetória caiu no buraco negro %zu %zu \n",j,i);
            }else{
                if((xf[1]>=6) && (xf[1]<=10)){
                    shadow[Nx*j + i] = 1.0;
                    /*double x_ap = cos(xf[3]);
                    double y_ap = sin(xf[3]);
                    if(x_ap>0){
                        
                        if(y_ap>0){
                            
                            shadow[Nx*j + i] = 1;
                        }else{
                            
                            shadow[Nx*j + i] = 2;
                        }
                    }else{
                         
                        if(y_ap>0){
                    
                            shadow[Nx*j + i] = 3;
                        }else{
                            
                            shadow[Nx*j + i] = 4;
                        }
                    }
                    //printf("Trajetória colidiu com disco %zu %zu \n",j,i);
                }else{
                    shadow[Nx*j + 1] = 0.0;
                    //printf("Trajetória escapou do buraco negro %zu %zu \n",j,i);
                }
                */
                }    
            }
        }
    }
    freeWL(ln);
    FILE *f = fopen("disco5.dat", "w");

    if (!f){
        printf("ERRO AO LER ARQUIVO\n");
        return 1;
    }
    for (size_t j = 0; j < Ny; j++) {
        for (size_t i = 0; i < Nx; i++) {
            fprintf(f, "%lf ", shadow[Nx*j + i]);
        }
        fprintf(f, "\n"); // quebra de linha para formato NxN
    }
    fclose(f);
    free(shadow);
    printf("Simulação finalizada");
}