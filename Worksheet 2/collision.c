#include "collision.h"
#include "LBDefinitions.h"

void computePostCollisionDistributions(
                                       double *currentCell, 
                                       const double * const tau, 
                                       const double *const feq
                                      ){
    
    int i;
    double tau_inv = 1.0 / *tau;
    
    for(i = 0; i < Q; i++){

        currentCell[i] = currentCell[i] - tau_inv * (currentCell[i] - feq[i]);

    }
    
}

void doCollision(
                 double *collideField, 
                 int *flagField,
                 const double * const tau,
                 int xlength
                ){
    
    double *currentCell;
    int x, y, z;
    int cellIdx;
    
    double density;
    double velocity[3];
    double feq[Q];
    
    int numGridPoints = xlength + 2;
    
    for(z = 1; z <= xlength; z++){
        
        for(y = 1; y <= xlength; y++){
            
            for(x = 1; x <= xlength; x++){
                
                //Determine current cell 
                cellIdx = Q * (z * numGridPoints * numGridPoints + y * numGridPoints + x);
                currentCell = collideField + cellIdx;
                
                //Compute post-collision distribution using density, velocity and equilibrium distribution 
                computeDensity(currentCell, &density);
                computeVelocity(currentCell, &density, velocity);
                computeFeq(&density, velocity, feq);
                computePostCollisionDistributions(currentCell, tau, feq);
   
            }

        }

    }
    
}


