#include "uvp.h"
#include <math.h>
#include <stdlib.h>
#include "helper.h"
#include "boundary_val.h"


/*Calculate the maximum time step allowed*/
void calculate_dt(
                  double Re,
                  double tau,
                  double *dt,
                  double dx,
                  double dy,
                  int il,
                  int ir,
                  int jb,
                  int jt,
                  double **U,
                  double **V
                  ){
    
    int i, j;
    double u_max = 0, v_max = 0;  
    double dt_min;
    
    /*Find the maximum velocities*/
    for (i = il; i <= ir; i++){
        
        for (j = jb; j <= jt; j++){
            
            if (fabs(U[i][j]) > u_max)
                
                u_max = fabs(U[i][j]);
            
            if (fabs(V[i][j]) > v_max)
                
                v_max = fabs(V[i][j]);
            
        }
        
    }
    
    /*Time step criteron*/
    dt_min = fmin(dx / u_max, dy / v_max);
    dt_min = fmin(0.5 * Re / (1 / (dx * dx) + 1 / (dy * dy)), dt_min);

    if (tau < 0){

	    *dt = *dt;

    }

    else{

    *dt = tau * dt_min;

    }

    
}


/*Calculate the F(n) and G(n)*/
void calculate_fg(
                  double Re,
                  double GX,
                  double GY,
                  double alpha,
                  double dt,
                  double dx,
                  double dy,
                  int il,
                  int ir,
                  int jb,
                  int jt,
                  int imax,
                  int jmax,
                  double **U,
                  double **V,
                  double **F,
                  double **G
                  ){
    
    int i, j;
    
        //Left wall BC
        if (il == 1){

            for(j = jb; j <= jt; j++){
            
            F[0][j] = U[0][j];
            
            }

        }

        //Right wall BC
        if (ir == imax){
        
            for(j = jb; j <= jt; j++){
    
                F[imax][j] = U[imax][j];

            }

        }

        /*Calculate F(n)*/
        for (i = max(1, il - 1); i <= min(ir, imax - 1); i++){
            
            for (j = jb; j <= jt; j++){
                
                F[i][j] = U[i][j] + dt * (
                                          + 1 / Re * ((U[i+1][j] - 2 * U[i][j] + U[i-1][j]) / (dx * dx) + (U[i][j+1] - 2 * U[i][j] + U[i][j-1]) / (dy * dy))
                                          - 1 / dx * (pow((U[i][j] + U[i+1][j]) / 2, 2) - pow((U[i-1][j] + U[i][j]) / 2, 2))
                                          - alpha / dx * (fabs(U[i][j] + U[i+1][j]) * (U[i][j] - U[i+1][j]) / 4 - fabs(U[i-1][j] + U[i][j]) * (U[i-1][j] - U[i][j]) / 4)
                                          - 1 / dy * ((V[i][j] + V[i+1][j]) * (U[i][j] + U[i][j+1]) / 4 - (V[i][j-1] + V[i+1][j-1]) * (U[i][j-1] + U[i][j]) / 4)
                                          - alpha / dy * (fabs(V[i][j] + V[i+1][j]) * (U[i][j] - U[i][j+1]) / 4 - fabs(V[i][j-1] + V[i+1][j-1]) * (U[i][j-1] - U[i][j]) / 4) + GX
                                          );
                
            }
            
        }
    
        //Lower wall BC
        if (jb == 1){

            for (i = il; i <= ir; i++){
        
                G[i][0] = V[i][0];
            
            }

        }
    
        //Upper wall BC
        if (jt == jmax){
        
            for (i = il; i <= ir; i++){
 
                G[i][jmax] = V[i][jmax];

            }

        }

        /*Calculate G(n)*/
        for (i = il; i <= ir; i++){
        
            for (j = max(1, jb - 1); j <= min(jt, jmax - 1); j++){
            
                G[i][j] = V[i][j] + dt * (
                                          + 1 / Re * ((V[i][j+1] - 2 * V[i][j] + V[i][j-1]) / (dy * dy) + (V[i+1][j] - 2 * V[i][j] + V[i-1][j]) / (dx * dx))
                                          - 1 / dy * (pow((V[i][j] + V[i][j+1]) / 2, 2) - pow((V[i][j-1] + V[i][j]) / 2, 2))
                                          - alpha / dy * (fabs(V[i][j] + V[i][j+1]) * (V[i][j] - V[i][j+1]) / 4 - fabs(V[i][j-1] + V[i][j]) * (V[i][j-1] - V[i][j]) / 4)
                                          - 1 / dx * ((U[i][j] + U[i][j+1]) * (V[i][j] + V[i+1][j]) / 4 - (U[i-1][j] + U[i-1][j+1]) * (V[i-1][j] + V[i][j]) / 4)
                                          - alpha / dx * (fabs(U[i][j] + U[i][j+1]) * (V[i][j] - V[i+1][j]) / 4 - fabs(U[i-1][j] + U[i-1][j+1]) * (V[i-1][j] - V[i][j]) / 4)
                                          + GY
                                          );
                
            }
            
        }
    
}


/*Calculate the RHS of the Pressure Poisson Equation*/
void calculate_rs(
                  double dt,
                  double dx,
                  double dy,
                  int il,
                  int ir,
                  int jb,
                  int jt,
                  int imax,
                  int jmax,
                  double **F,
                  double **G,
                  double **RS
                  ){
    
    int i, j;
    
    for (i = il; i <= ir; i++){
        
        for (j = jb; j <= jt; j++){
            
            RS[i][j] = 1 / dt * ((F[i][j] - F[i-1][j]) / dx + (G[i][j] - G[i][j-1]) / dy);
            
        }
    }
    
}


/*Update the velocities*/
void calculate_uv(
                  double dt,
                  double dx,
                  double dy,
                  int il,
                  int ir,
                  int jb,
                  int jt,
                  int imax,
                  int jmax,
                  double **U,
                  double **V,
                  double **F,
                  double **G,
                  double **P
                  ){
    
    int i, j;
    
        /*Update the U velocities*/
        for (i = max(1, il - 1); i <= min(ir, imax - 1); i++){
            
            for (j = jb; j <= jt; j++){
                
                U[i][j] = F[i][j] - dt * (P[i+1][j] - P[i][j]) / dx;
                
            }
            
        }
    
        /*Update the V velocities*/
        for (i = il; i <= ir; i++){
            
            for (j = max(1, jb - 1); j <= min(jt, jmax - 1); j++){
                
                V[i][j] = G[i][j] - dt * (P[i][j+1] - P[i][j]) / dy;
                
            }
        }
    
}


