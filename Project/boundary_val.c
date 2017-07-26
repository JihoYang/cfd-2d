#include "boundary_val.h"
#include <stdio.h>
#include "init.h" 
#include "helper.h"

/*Set the Boundary Conditions (velocities)*/
void boundaryvalues(
                     int imax,
                     int jmax,
                     double **U,
                     double **V,
                     int wl,
                     int wr,
                     int wt,
                     int wb,
                     int **Flag
                   ){
    
    // no need to extract and external pgm file for this condition to be met
    
    int i,j;
    
    //Left boundary
    switch(wl){
            
        //No Slip
        case 1:
            
            for(j = 0; j <= jmax; j++){
                
                U[0][j] = 0;
                V[0][j] = - V[1][j];
                
            }
            
            break;
            
        //Free Slip
        case 2:
            
            for(j = 0; j <= jmax; j++){
                
                U[0][j] = 0;
                V[0][j] = V[1][j];
                
            }
            
            break;
            
        //Outflow
        case 3:
            
            for(j = 0; j <= jmax; j++){
                
                U[0][j] = U[1][j];
                V[0][j] = V[1][j];
                
            }
            
            break;
            
    }
    
    
    //Right boundary
    switch(wr){
            
        //No Slip
        case 1:
            
            for(j = 0; j <= jmax; j++){
                
                U[imax][j] = 0;
                V[imax+1][j] = -V[imax][j];
                
            }
            
            break;
            
        //Free Slip
        case 2:
            
            for(j = 0; j <= jmax; j++){
                
                U[imax][j] = 0;
                V[imax+1][j] = V[imax][j];
                
            }
            
            break;
            
        //Outflow
        case 3:
            
            for(j = 0; j <= jmax; j++){
                
                U[imax][j] = U[imax-1][j];
                V[imax+1][j] = V[imax][j];
                
            }
            
            break;
            
    }
    
    
    //Top boundary
    switch(wt){
            
        //No Slip
        case 1:
            
            for(i = 0; i <= imax; i++){
                
                U[i][jmax+1] = -U[i][jmax];
                V[i][jmax] = 0;
                
            }
            
            break;
            
        //Free Slip
        case 2:
            
            for(i = 0; i <= imax; i++){
                
                U[i][jmax+1] = U[i][jmax];
                V[i][jmax] = 0;
                
            }
            
            break;
            
        //outflow
        case 3:
            
            for(i = 0; i <= imax; j++){
                
                U[i][jmax+1] = U[i][jmax];
                V[i][jmax] = V[i][jmax-1];
                
            }
            
            break;
            
        //Driven cavity flow
        case 4:
            
            for( i=0; i<= imax; i++){
                
                U[i][jmax+1] = 2* 1.0  - U[i][jmax];
                V[i][jmax] =0;
                
            }
            
            break;
            
            
    }
    
    //Bottom boundary
    switch(wb){
            
            
        //No Slip
        case 1:
            
            for(i = 0; i <= imax; i++){
                
                U[i][0] = -U[i][1];
                V[i][0] = 0;
                
            } 
            
            break;
            
            
        //Free Slip
        case 2:
            
            for(i = 0; i <= imax; j++){
                
                U[i][0] = U[i][1];
                V[i][0] = 0;
                
            }
            
            break;
            
            
        //Outflow
        case 3:
            
            for(i = 0; i <= imax; j++){
                
                U[i][0] = U[i][1];
                V[i][0] = V[i][1];
                
            }
            
            break;
            
    }
    
    //Arbitrary geometry (obstacles)
    for(i = 1; i <= imax; i++){
        
        for(j = 1; j <= jmax; j++){
            
            switch(Flag[i][j]){
                    
                case B_N:
                    
                    V[i][j] = 0;
                    U[i-1][j] = -U[i-1][j+1];
                    U[i][j] = -U[i][j+1];
                    
                    break;
                    
                case B_S:
                    
                    V[i][j-1]=0;
                    U[i-1][j] = -U[i-1][j-1];
                    U[i][j] = -U[i][j-1];
                    
                    break;
                    
                case B_O:
                    
                    U[i][j] = 0;
                    V[i][j] = -V[i+1][j];
                    V[i][j-1] = -V[i+1][j-1];
                    
                    break;
                    
                case B_W:
                    
                    U[i-1][j] = 0;
                    V[i][j] = -V[i-1][j];
                    V[i][j-1] = -V[i-1][j-1];
                    
                    break;
                    
                case B_NW:
                    
                    V[i][j] = 0;
                    V[i][j-1] = -V[i-1][j-1];
                    U[i-1][j] = 0;
                    U[i][j] = -U[i][j+1];
                    
                    break;
                    
                case B_NO:
                    
                    U[i][j] = 0;
                    V[i][j] = 0;
                    U[i-1][j] = -U[i-1][j+1];
                    V[i][j-1] = -V[i+1][j-1];
                    
                    break;
                    
                case B_SW:
                    
                    U[i-1][j] = 0;
                    V[i][j-1] = 0;
                    U[i][j] = -U[i][j-1];
                    V[i][j] = -V[i-1][j];
                    
                    break;
                    
                case B_SO:
                    
                    U[i][j] = 0;
                    V[i][j-1] = 0;
                    V[i][j] = -V[i+1][j];
                    U[i-1][j] = -U[i-1][j-1];
                    
                    break;
                    
            }
            
        }
        
    }
    
}

