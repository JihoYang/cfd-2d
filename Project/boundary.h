#ifndef __BOUNDARY__H__
#define __BOUNDARY__H__
 
void spec_boundary_val(
                        char *problem, // Defines the problem type in a character or character string
                        int imax,      // Grid length along x direction
                        int jmax,      // Grid length along y direction 
                        double **U,    // U maxtrix, x component of velocity 
                        double **V,    // V matrix, y component of velocity
                        double dP,     // Pressure difference
                        double Re,     // Reynodls number
                        double dx,     // grid size - x direction
                        double dy,     // grid size - y direction
                        int **Flag     // Flagfield
                      );


#endif
