#include "computeCellValues.h"
#include "LBDefinitions.h"
#include <math.h>

void computeDensity(
                    const double *const currentCell, 
                    double *density
                   ){

	int i;
	
	*density = 0;
	
	for(i = 0; i < Q; i++){

        *density += currentCell[i];

	}

}

void computeVelocity(
                     const double * const currentCell, 
                     const double * const density,  
                     double *velocity
                    ){

	int i;
	
	velocity[0] = 0;
	velocity[1] = 0;
	velocity[2] = 0;
	
	for(i = 0; i < Q; i++){

        velocity[0] += currentCell[i] * LATTICEVELOCITIES[i][0];
		velocity[1] += currentCell[i] * LATTICEVELOCITIES[i][1];
		velocity[2] += currentCell[i] * LATTICEVELOCITIES[i][2];

	}

	velocity[0] = velocity[0] / *density;
	velocity[1] = velocity[1] / *density;
	velocity[2] = velocity[2] / *density;

}

void computeFeq(
                const double * const density, 
                const double * const velocity, 
                double *feq
               ){

	int i;
	
	double u_dot_u = velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2];
	double div_C_S_2 = 1.0 / (C_S * C_S);
	double div_C_S_4 = div_C_S_2 * div_C_S_2;
	
	for(i = 0; i < Q; i++){

		double c_dot_u = LATTICEVELOCITIES[i][0] * velocity[0] + LATTICEVELOCITIES[i][1] * velocity[1] + LATTICEVELOCITIES[i][2] * velocity[2];
		double T1 = c_dot_u * div_C_S_2;
		double T2 = 0.5 * c_dot_u * c_dot_u * div_C_S_4;
		double T3 = 0.5 * u_dot_u * div_C_S_2;

		feq[i] = LATTICEWEIGHTS[i] * (*density) * (1 + T1 + T2 - T3);

	}

}

