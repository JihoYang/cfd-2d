#include "initLB.h"
#include <math.h>
#include "LBDefinitions.h"

int readParameters(
                   const char *szFileName, 
                   int *xlength, 
                   double *velocityWall, 
                   double *Re, 
                   int *timesteps, 
                   int *timestepsPerPlotting, 
                   int argc, 
                   char *argv[]
                  ){
    
    READ_INT(szFileName, *xlength);
    READ_DOUBLE(szFileName, *velocityWall);
    READ_DOUBLE(szFileName, *Re);
    READ_INT(szFileName, *timesteps);
    READ_INT(szFileName, *timestepsPerPlotting);
    
    return 0;

}


void initialiseFields(
                      double *collideField, 
                      double *streamField, 
                      int *flagField, 
                      int xlength
                     ){
 
	int x,y,z,i;
	int cellIdx;
	int numGridPoints = xlength + 2;

	//Initialize flagField to indicate type of cell: FLUID=0, NO SLIP=1 and MOVING WALL=2
	for(z = 0; z <= xlength + 1; z++){

	    for(y = 0; y <= xlength + 1; y++){

		    for(x = 0; x <= xlength + 1; x++){
                    
                cellIdx = z * numGridPoints * numGridPoints + y * numGridPoints + x ;

                if(z == xlength + 1){

                    flagField[cellIdx] = MOVING_WALL;

                }

                else if(x==0 || y==0 || z==0 || x==xlength+1 || y==xlength+1){

                    flagField[cellIdx] = NO_SLIP;

                }

                else{

                    flagField[cellIdx] = FLUID;

                }

            }
            
        }

    }

	//Initialize matrices of collision and stream: f(x,t=0) = omega_i
	for(z = 0; z <= xlength + 1; z++){

	    for(y = 0; y <= xlength + 1; y++){

	        for(x = 0; x <= xlength + 1; x++){
                
			    for(i = 0; i < Q; i++){
                    
                    cellIdx = Q * (z * numGridPoints * numGridPoints + y * numGridPoints + x) + i;
                    
					collideField[cellIdx] = LATTICEWEIGHTS[i];
					streamField[cellIdx] = LATTICEWEIGHTS[i];
                    
                }

            }

        }

    }

}

