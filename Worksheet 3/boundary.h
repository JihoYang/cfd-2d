#ifndef __BOUNDARY__H__
#define __BOUNDARY__H__
 
void spec_boundary_val(
                        char *problem, /* defines the problem type in a character or character string */
                        int imax,      /* grid length along x direction */
                        int jmax,      /* grid length along y direction */
                        double **U,    /* U maxtrix, x component of velocity */
                        double **V,     // V matrix, y component of velocity
                        double dP,
                        double Re,
                        double dx,
                        double dy
                       );


#endif
