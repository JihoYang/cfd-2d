#ifndef _VISUALLB_H_
#define _VISUALLB_H_

//Writes the density and velocity Field (derived from the distributions in collideField) to a file determined by 'filename' and timestep 't'
void writeVtkOutput(
                    const double * const collideField, 
                    const int * const flagField, 
                    const char * filename, 
                    unsigned int tau, 
                    int xlength
);

void write_vtkHeader( 
                     FILE *fp, 
                     int xlength 
);

void write_vtkPointCoordinates( 
                               FILE *fp, 
                               int numGridPoints
);

#endif

