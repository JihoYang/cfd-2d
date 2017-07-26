#include "streaming.h"
#include "LBDefinitions.h"
#include <stdlib.h>
#include <stdio.h>

void doStreaming(
                 double *collideField, 
                 double *streamField,
                 int *flagField,
                 int xlength
                ){
    
    //Index for i-th direction 
    int i;
    
    //Coordinates of cell to be streamed to 
    int x, y, z;
    
    //Coordinates of neighbour cells 
    int xn, yn, zn;
    
    //Pointer and array indices for current and neighbor cell 
    double *currentCell, *neighborCell;
    int cell_index, neighborCell_index;
    
    //Total number of grid points per side 
    int numGridPoints = xlength + 2;
    
    //Loop through all fluid cells
    for(z = 1; z <= xlength; z++){
        
        for(y = 1; y <= xlength; y++){
            
            for(x = 1; x <= xlength; x++){
                
                //Determine current cell at (x, y, z) 
                cell_index = Q * (z * numGridPoints * numGridPoints + y * numGridPoints + x);
                currentCell = streamField + cell_index;
                
                for(i = 0; i < Q; i++){
                
                    xn = x - LATTICEVELOCITIES[i][0];
                    yn = y - LATTICEVELOCITIES[i][1];
                    zn = z - LATTICEVELOCITIES[i][2];

                    neighborCell_index = Q * (zn * numGridPoints * numGridPoints + yn * numGridPoints + xn);
                    neighborCell = collideField + neighborCell_index;
           
                    currentCell[i] = neighborCell[i];
                    
                }
                
            }	

        }

    }

}

