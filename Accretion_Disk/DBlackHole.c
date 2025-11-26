#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include "DBlackHole.h"

struct line{
    double x0[6]; // Condições iniciais
    double x[6];  // Condições atuais 
    double m; double q; double E; double lz; double Qcarter; double M; double a; double Q;
    double dt;
    double N;
};

Linha_de_mundo* create_Wline(size_t N, double dt,double* x0, double m, double q, double E, double lz, double Qcarter, double M, double Q, double a){

    Linha_de_mundo* ln = (Linha_de_mundo*) malloc(sizeof(Linha_de_mundo));
    for( size_t i = 0; i<6 ; i++){
        ln->x0[i] = x0[i];
        ln->x[i] = x0[i];
    }
    ln->m = m;
    ln->q = q;
    ln->E = E; 
    ln->lz = lz; 
    ln->Qcarter = Qcarter;
    ln->M = M; ln->Q = Q; ln->a = a;

    ln->N = N;
    ln->dt = dt;
    if(!isnan(Qcarter)){
        // A escolha de sinal é necessária em ambos os casos
            double cosseno = cos(x0[2]);
            double seno = sin(x0[2]);
            double delta = x0[1]*x0[1] + ln->a*ln->a + ln->Q*ln->Q - 2*ln->M*x0[1];
            double util = ((ln->a*ln->a + x0[1]*x0[1])*ln->E - ln->a*ln->lz - ln->q*ln->Q*x0[1] );
            //pivo = (-Qcarter - (lz-a*E)**2 - m*m*r*r + (util*util/delta) )/delta
            double pivo = (-Qcarter - (lz-a*E)*(lz-a*E) - m*m*x0[1]*x0[1] + (util*util/delta) )/delta;
            if((pivo<0) && (fabs(pivo)>pow(10,-30))){
                printf("ERRO: Trajetória proibida R(r)<0 ");
                printf("Condição Inicial será p_r = 0 \n");
                ln->x0[4] = 0.0;
                ln->x[4] = 0.0; 
            }else{
                ln->x0[4] = (x0[4]/fabs(x0[4])) * sqrt(pivo) ;
                ln->x[4] = ln->x0[4];
            }
            pivo = Qcarter + (lz-a*E)*(lz-a*E) - (((lz-a*E*seno*seno)/seno)*((lz-a*E*seno*seno)/seno)) - m*m*a*a*cosseno*cosseno;
            if((pivo<0) && (fabs(pivo)>pow(10,-30)) ){
                 printf("ERRO: Trajetória proibida Theta(theta)<0 ");
                printf("Condição Inicial será p_theta = 0\n");
                ln->x0[5] = 0.0;
                ln->x[5] = 0.0; 
            }else{
                ln->x0[5] = (x0[5]/fabs(x0[5])) * sqrt(pivo) ;
                ln->x[5] = ln->x0[5];
            }
    }

    return ln;
}

Linha_de_mundo* reuse_Wline(Linha_de_mundo* ln,size_t N, double dt,double* x0, double m, double q, double E, double lz, double Qcarter, double M, double Q, double a){

    for( size_t i = 0; i<6 ; i++){
        ln->x0[i] = x0[i];
        ln->x[i] = x0[i];
    }
    ln->m = m;
    ln->q = q;
    ln->E = E; 
    ln->lz = lz; 
    ln->Qcarter = Qcarter;
    ln->M = M; ln->Q = Q; ln->a = a;

    ln->N = N;
    ln->dt = dt;
    if(!isnan(Qcarter)){
        // A escolha de sinal é necessária em ambos os casos
            double cosseno = cos(x0[2]);
            double seno = sin(x0[2]);
            double delta = x0[1]*x0[1] + ln->a*ln->a + ln->Q*ln->Q - 2*ln->M*x0[1];
            double util = ((ln->a*ln->a + x0[1]*x0[1])*ln->E - ln->a*ln->lz - ln->q*ln->Q*x0[1] );
            //pivo = (-Qcarter - (lz-a*E)**2 - m*m*r*r + (util*util/delta) )/delta
            double pivo = (-Qcarter - (lz-a*E)*(lz-a*E) - m*m*x0[1]*x0[1] + (util*util/delta) )/delta;
            if((pivo<0) && (fabs(pivo)>pow(10,-30))){
                printf("ERRO: Trajetória proibida R(r)<0 ");
                printf("Condição Inicial será p_r = 0 \n");
                ln->x0[4] = 0.0;
                ln->x[4] = 0.0; 
            }else{
                ln->x0[4] = (x0[4]/fabs(x0[4])) * sqrt(pivo) ;
                ln->x[4] = ln->x0[4];
            }
            pivo = Qcarter + (lz-a*E)*(lz-a*E) - (((lz-a*E*seno*seno)/seno)*((lz-a*E*seno*seno)/seno)) - m*m*a*a*cosseno*cosseno;
            if((pivo<0) && (fabs(pivo)>pow(10,-30)) ){
                 printf("ERRO: Trajetória proibida Theta(theta)<0 ");
                printf("Condição Inicial será p_theta = 0\n");
                ln->x0[5] = 0.0;
                ln->x[5] = 0.0; 
            }else{
                ln->x0[5] = (x0[5]/fabs(x0[5])) * sqrt(pivo) ;
                ln->x[5] = ln->x0[5];
            }
    }

    return ln;
}







void dudt(Linha_de_mundo* ln, double* x, double* k){
    
        //t, r , theta , phi , pr , ptheta 
        //E, lz,m,q,M,Q,a = const
        // Expressões algébricas úteis
        double cosseno = cos(x[2]);
        double seno = sin(x[2]);
        double delta = x[1]*x[1] + ln->a*ln->a + ln->Q*ln->Q - 2*ln->M*x[1];
        double sigma = x[1]*x[1] + ln->a*ln->a*(cosseno*cosseno);
        double util = ((ln->a*ln->a + x[1]*x[1])*ln->E - ln->a*ln->lz - ln->q*ln->Q*x[1] );
        
        if(delta<1e-5){
            // Partícula possivelmente já caiu dentro do buraco negro, dinâmica explodindo, impede dinâmica de continuar.
            k[0] = -10;
            return;
        }
        
        k[0] = (1/sigma)*(-ln->a*(ln->a*ln->E*seno*seno - ln->lz)   +  ((x[1]*x[1] + ln->a*ln->a)*(util)/delta)    ) ;
    
        
        k[1] = (delta/sigma) * x[4] ;

         
        k[2] = (1/sigma) * x[5];
        
        
        k[3] = (1/sigma)*( -(ln->a*ln->E - ln->lz/(seno*seno)) + (ln->a/delta)*util  );

        k[4] = (1/(sigma))*(-((util*util*(-ln->M+x[1]))/(delta*delta)) + (util*(2*x[1]*ln->E - ln->q*ln->Q)/(delta)) - ((-ln->M + x[1])*x[4]*x[4]) - (ln->m*ln->m*x[1])    );
        
        
        k[5] =  (1/sigma)*( cosseno*seno*(ln->a*ln->a*(-ln->E*ln->E + ln->m*ln->m) + (ln->lz*ln->lz/(seno*seno)) ) + ((cosseno*cosseno*cosseno)*ln->lz*ln->lz/(seno*seno*seno)) );
        
}

double* evolve_RK(Linha_de_mundo* ln){
    for(size_t j = 0; j<ln->N; j++){
        double k1[6];
        double k2[6];
        double k3[6];
        double k4[6];
        double temp[6];
        dudt(ln,ln->x,k1); // Estimativa para K1
    
        for(size_t i = 0; i<6; i++){temp[i] = ln->x[i] + (k1[i]*ln->dt/2);}
        dudt(ln,temp,k2); // Estimativa para K2

        for(size_t i = 0; i<6; i++){temp[i] = ln->x[i] + (k2[i]*ln->dt/2);}
        dudt(ln,temp,k3); // Estimativa para K3

        for(size_t i = 0; i<6; i++){temp[i] = ln->x[i] + (k3[i]*ln->dt);}
        dudt(ln,temp,k4); // Estimativa para K4

        // Atualiza posição e trajetórias das órbitas 
        for(size_t i = 0; i<6; i++){ln->x[i] = ln->x[i] + (ln->dt/6)*(k1[i]+ 2*k2[i] + 2*k3[i] + k4[i]);}
        if((ln->x[1]<1.0001) || (k1[1]<0)){
            break;
        }
    }
    return ln->x;
    // Fim da implementação da trajetória    
}

double* evolve_RK_inverse(Linha_de_mundo* ln){
    // Raio do horizonte de eventos
    double r_bh = ln->M + sqrt(ln->M*ln->M - ln->a*ln->a -ln->Q*ln->Q);
    double k1[6];
    double k2[6];
    double k3[6];
    double k4[6];
    double temp[6]; 
    for(size_t j = 0; j<ln->N; j++){
        // Testa se a partícula foi capturada
        if(ln->x[1]<=r_bh*1.0001){
            // Partícula capturada
            //printf("Partícula capturada pois chegou muito próximo do buraco negro");
            break;
        }
        dudt(ln,ln->x,k1); // Estimativa para K1
    
        for(size_t i = 0; i<6; i++){temp[i] = ln->x[i] - (k1[i]*ln->dt/2);}
        dudt(ln,temp,k2); // Estimativa para K2

        for(size_t i = 0; i<6; i++){temp[i] = ln->x[i] - (k2[i]*ln->dt/2);}
        dudt(ln,temp,k3); // Estimativa para K3

        for(size_t i = 0; i<6; i++){temp[i] = ln->x[i] - (k3[i]*ln->dt);}
        dudt(ln,temp,k4); // Estimativa para K4

        // Atualiza posição e trajetórias das órbitas 
        for(size_t i = 0; i<6; i++){temp[i] = ln->x[i] - (ln->dt/6)*(k1[i]+ 2*k2[i] + 2*k3[i] + k4[i]);}

        double rho_prev = ln->x[1]*sin(ln->x[2]);
        double z_prev = ln->x[1]*cos(ln->x[2]);

        double rho = temp[1]*sin(temp[2]);
        double z = temp[1]*cos(temp[2]);

        if(z_prev*z<0){
            // Partícula atravessou o equador
            double rho0 = rho_prev - (rho-rho_prev)*(z_prev)/(z-z_prev);
            if(((rho0>=6.0) && (rho0<=10.0))){
                // A Luz veio do disco de acreção
                ln->x[1] = rho0;
                return ln->x;
            }
        }

        for(size_t i = 0; i<6; i++){ln->x[i] = temp[i];}
        
    
    }
    return ln->x;
    // Fim da implementação da trajetória    
}


size_t freeWL(Linha_de_mundo* ln){
    free(ln);
    return 0;
}

// Retorna posição da Linha_de_mundo
double* get_x(Linha_de_mundo* ln){
    return ln->x;
}

// Retorna posição original da Linha_de_mundo
double* get_x0(Linha_de_mundo* ln){
    return ln->x0;
}