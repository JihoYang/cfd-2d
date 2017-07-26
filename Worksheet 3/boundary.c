#include "boundary.h"
#include "init.h"
#include <stdio.h> 
#include <string.h>


void spec_boundary_val (char *problem, int imax, int jmax, double **U, double **V, double dP, double Re, double dx, double dy){
            
            double U_max = 1.0;
            int i, j;
            
            if(dP == 0){

                if(strcmp(problem, "karman_vortex_street") == 0){

                    for(j = 1; j <= jmax; j++){

                        U[0][j] = 1;
                        V[0][j] = -V[1][j];

                    }

                }

                else if(strcmp(problem, "plane_shear_flow") == 0){
    
                    for(j = 1 ; j <= jmax; j++) {

                        //U[0][j] = U_max * (1 - j*j / (jmax * jmax));
                        U[0][j] = 1;
                        V[0][j] = -V[1][j];     

                    }

                }
        
                else if(strcmp(problem, "flow_over_step") == 0){

                    for(j  = jmax/2 + 1; j <= jmax; j++){

                        U[0][j] = 1;
                        V[0][j] = -V[1][j];

                    }
    
                }

                else if(strcmp(problem, "driven_cavity") == 0){

                    for(i = 1; i <= imax; i++){

                        U[i][jmax + 1] = 2 * U_max - U[i][jmax];
                        V[i][jmax] = 0;

                    }

                } 

            }
    
            if(dP != 0){
        
                if(strcmp(problem, "plane_shear_flow") == 0){
            
                    for(j = 1; j <= jmax; j++) {
                
                        U[0][j] = -0.5 * Re * dP / dx * dy * (dy - jmax);
                        V[0][j] = 0; 

                    }
            
                }
                
            }


}
            
            
