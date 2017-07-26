#include "boundary.h"
#include "init.h"
#include <stdio.h> 
#include <string.h>

void spec_boundary_val (char *problem, int imax, int jmax, double **U, double **V, double dP, double Re, double dx, double dy, int **Flag){
    
    int j;
    //double U_max = 1.0;
    
    //Inflow conditions (Dirichlet velocity BC)
    if(strcmp(problem, "layout1") == 0) {

        for(j = 0; j <= jmax; j++) {

            U[0][j] = 1;
            V[0][j] = 0;

        }

    } 

    else if(strcmp(problem, "layout2") == 0){

        for(j = 0; j <= jmax; j++){

            U[0][j] = 1;
            V[0][j] = 0;

        }

    }

    else if(strcmp(problem, "layout3") == 0){

        for(j = 0; j <= jmax; j++){

            U[0][j] = 1;
            V[0][j] = 0;

        }
  }
  
    else if(strcmp(problem, "layout4") == 0){
      
	for(j = 0; j <= jmax; j++){

            U[0][j] = 1;
            V[0][j] = 0;
       }

  }	
	
   else if(strcmp(problem, "layout5") == 0){

        for(j = 0; j <= jmax; j++){

            U[0][j] = 1;
            V[0][j] = 0;

      }

  }


}


            
