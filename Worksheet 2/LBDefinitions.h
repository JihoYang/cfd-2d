#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_

#include <math.h>

//Dimension and particle direction components
#define Q 			19
#define D 			3

//Lattice Weights
#define W_12_36		12.0/36.0
#define W_2_36		2.0/36.0
#define W_1_36		1.0/36.0

//Flags
#define FLUID			0
#define NO_SLIP		    1
#define MOVING_WALL	    2

//Speed of sound (1 / sqrt(3))
static const double C_S = 0.57735026919;

//Lattice velocities
static const int LATTICEVELOCITIES[19][3] = {
	
    { 0, -1, -1},
    {-1,  0, -1},
    { 0,  0, -1},
    { 1,  0, -1},
    { 0,  1, -1},
    
    {-1, -1,  0},
    { 0, -1,  0},
    { 1, -1,  0},
    {-1,  0,  0},
    { 0,  0,  0},
    
    { 1,  0,  0},
    {-1,  1,  0},
    { 0,  1,  0},
    { 1,  1,  0},
    { 0, -1,  1},
    
    {-1,  0,  1},
    { 0,  0,  1},
    { 1,  0,  1},
    { 0,  1,  1}
    
};


//Lattice weights
static const double LATTICEWEIGHTS[19] = {

    W_1_36,
    W_1_36,
    W_2_36,
    W_1_36,
    W_1_36,
    
    W_1_36,
    W_2_36,
    W_1_36,
    W_2_36,
    W_12_36,
    
    W_2_36,
    W_1_36,
    W_2_36,
    W_1_36,
    W_1_36,
    
    W_1_36,
    W_2_36,
    W_1_36,
    W_1_36		
    
};

#endif

