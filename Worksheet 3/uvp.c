#include "uvp.h"
#include <math.h>
#include <stdlib.h>
#include "helper.h"
#include "init.h"


/*Calculate the maximum time step allowed*/
void calculate_dt(
                  double Re,
                  double tau,
                  double *dt,
                  double dx,
                  double dy,
                  int imax,
                  int jmax,
                  double **U,
                  double **V
                  ){
    
    int i, j;
    double u_max = U[1][1], v_max = V[1][1];
    double dt_min;
    
    /*Find the maximum velocities*/
    for (i = 1; i <= imax; i++){
        
        for (j = 1; j <= jmax; j++){
            
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
                  int imax,
                  int jmax,
                  double **U,
                  double **V,
                  double **F,
                  double **G,
                  int** Flag
                  ){
    
    int i, j;

    //Boundary condition F and G at left wall (taken from WS1)
    for(j = 1; j <= jmax; j++){
            
        F[0][j] = U[0][j];
        F[imax][j] = U[imax][j];
            
    }    

    for (i = 1; i <= imax; i++){
         
        G[i][0] = V[i][0];                                                                                                                                                   
        G[i][jmax] = V[i][jmax];
 
    }

    //Boundary condition for inner boundaries
    for(i = 1; i <= imax; i++){
    
        for( j  = 1; j <= jmax;j++){
    
            if(Flag[i][j] == B_O){

            F[i][j] = U[i][j];

            }
         
            else if(Flag[i][j] == B_SO){     

                F[i][j] = U[i][j];      
                G[i][j-1] = V[i][j-1];

            }

            else if(Flag[i][j] == B_S){

                G[i][j-1] = V[i][j-1];
    
            }
         
            else if(Flag[i][j] == B_SW){
    
                G[i][j-1] = V[i][j-1];
                F[i-1][j] = U[i-1][j];

            } 

            else if(Flag[i][j] == B_NO){

                G[i][j] = V[i][j];
                F[i][j] = U[i][j];

            }

            else if(Flag[i][j] == B_W){
          
            F[i-1][j] = U[i-1][j];

            }
        
            else if(Flag[i][j] == B_NW){

            G[i][j] = V[i][j];
            F[i-1][j] = -U[i-1][j];
    
            }
  
            else if(Flag[i][j] == B_N){

            G[i][j] = V[i][j];

            }

        }

    }

    /*Calculate F(n)*/
    for (i = 1; i <= imax - 1; i++){
            
        for (j = 1; j <= jmax; j++){
                
            F[i][j] = U[i][j] + dt * (
                                      + 1 / Re * ((U[i+1][j] - 2 * U[i][j] + U[i-1][j]) / (dx * dx) + (U[i][j+1] - 2 * U[i][j] + U[i][j-1]) / (dy * dy))
                                      - 1 / dx * (pow((U[i][j] + U[i+1][j]) / 2, 2) - pow((U[i-1][j] + U[i][j]) / 2, 2))
                                      - alpha / dx * (fabs(U[i][j] + U[i+1][j]) * (U[i][j] - U[i+1][j]) / 4 - fabs(U[i-1][j] + U[i][j]) * (U[i-1][j] - U[i][j]) / 4)
                                      - 1 / dy * ((V[i][j] + V[i+1][j]) * (U[i][j] + U[i][j+1]) / 4 - (V[i][j-1] + V[i+1][j-1]) * (U[i][j-1] + U[i][j]) / 4)
                                      - alpha / dy * (fabs(V[i][j] + V[i+1][j]) * (U[i][j] - U[i][j+1]) / 4 - fabs(V[i][j-1] + V[i+1][j-1]) * (U[i][j-1] - U[i][j]) / 4) + GX
                                      );
                
        }
          
    }
    
    /*Calculate G(n)*/
    for (i = 1; i <= imax; i++){
        
        for (j = 1; j <= jmax - 1; j++){
            
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
                  int imax,
                  int jmax,
                  double **F,
                  double **G,
                  double **RS
                  ) {
    
    int i, j;
    
    for (i = 1; i <= imax; i++){
        
        for (j = 1; j <= jmax; j++){
            
            RS[i][j] = 1 / dt * ((F[i][j] - F[i-1][j]) / dx + (G[i][j] - G[i][j-1]) / dy);
            
        }
    }
    
}


/*Update the velocities*/
void calculate_uv(
                  double dt,
                  double dx,
                  double dy,
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
        for (i = 1; i <= imax - 1; i++){
            
            for (j = 1; j <= jmax; j++){
                
                U[i][j] = F[i][j] - dt * (P[i+1][j] - P[i][j]) / dx;
                
            }
            
        }
    
        /*Update the V velocities*/
        for (i = 1; i <= imax; i++){
            
            for (j = 1; j <= jmax - 1; j++){
                
                V[i][j] = G[i][j] - dt * (P[i][j+1] - P[i][j]) / dy;
                
            }
        }
    
}


