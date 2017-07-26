#include "boundary_val.h"
#include <stdio.h>
#include "helper.h"

/*Set the Boundary Conditions (velocities)*/
void boundaryvalues(
                    int il,
                    int ir,
                    int jb,
                    int jt,
                    int imax,
                    int jmax,
                    double **U,
                    double **V
                    ){
    
    int i, j;
    double U_mv_wall = 1.0;
    
    /*BC @ left wall*/
    if(il == 1){

        for(j = jb; j <= jt; j++){
            
            U[0][j] = 0;
            V[0][j] = -V[1][j];
        
        }

    }
    
    /*BC @ right wall*/
    if(ir == imax){

        for(j = jb; j <= jt; j++) {
        
            U[imax][j] = 0;
            V[imax+1][j] = -V[imax][j];
        
        }

    }
    
    /*BC @ floor*/
    if (jb == 1){

        for(i = il; i <= ir; i++) {
        
            U[i][0] = -U[i][1];
            V[i][0] = 0;
        
        }

    }
   
    /*BC @ moving wall*/
    if (jt == jmax){ 

    for(i = il; i <= ir; i++) {
        
        U[i][jmax+1] = 2 * U_mv_wall - U[i][jmax];
        V[i][jmax] = 0;
        
    }

    }
    
}

