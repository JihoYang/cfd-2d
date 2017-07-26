#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

//Handles the boundaries in our simulation setup
void treatBoundary(
                   double *collideField, 
                   int* flagField, 
                   const double * const wallVelocity,
                   int xlength
);

#endif

