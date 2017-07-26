#include "helper.h"
#include "init.h"
#include <mpi.h>

int read_parameters( const char *szFileName,       /* name of the file */
                    double *Re,                    /* reynolds number   */
                    double *UI,                    /* velocity x-direction */
                    double *VI,                    /* velocity y-direction */
                    double *PI,                    /* pressure */
                    double *GX,                    /* gravitation x-direction */
                    double *GY,                    /* gravitation y-direction */
                    double *t_end,                 /* end time */
                    double *xlength,               /* length of the domain x-dir.*/
                    double *ylength,               /* length of the domain y-dir.*/
                    double *dt,                    /* time step */
                    double *dx,                    /* length of a cell x-dir. */
                    double *dy,                    /* length of a cell y-dir. */
                    int  *imax,                    /* number of cells x-direction*/
                    int  *jmax,                    /* number of cells y-direction*/
                    double *alpha,                 /* uppwind differencing factor*/
                    double *omg,                   /* relaxation factor */
                    double *tau,                   /* safety factor for time step*/
                    int  *itermax,                 /* max. number of iterations  */
                                                   /* for pressure per time step */
                    double *eps,                   /* accuracy bound for pressure*/
                    double *dt_value,              /* time for output */
                    int *iproc,                    // number of subdomain cells x-direction
                    int *jproc,                    // number of subdomain cells y-direction
                    int myrank                     // current rank. here it only acts as switch for printing texts
                    )
{
    
    READ_DOUBLE( szFileName, *xlength, myrank );
    READ_DOUBLE( szFileName, *ylength, myrank );
    
    READ_DOUBLE( szFileName, *Re    , myrank);
    READ_DOUBLE( szFileName, *t_end , myrank);
    READ_DOUBLE( szFileName, *dt    , myrank);
    
    READ_INT   ( szFileName, *imax , myrank);
    READ_INT   ( szFileName, *jmax , myrank);
    
    READ_DOUBLE( szFileName, *omg   , myrank);
    READ_DOUBLE( szFileName, *eps   , myrank);
    READ_DOUBLE( szFileName, *tau   , myrank);
    READ_DOUBLE( szFileName, *alpha , myrank);
    
    READ_INT   ( szFileName, *itermax , myrank);
    READ_DOUBLE( szFileName, *dt_value, myrank);
    
    READ_DOUBLE( szFileName, *UI , myrank);
    READ_DOUBLE( szFileName, *VI , myrank);
    READ_DOUBLE( szFileName, *GX , myrank);
    READ_DOUBLE( szFileName, *GY , myrank);
    READ_DOUBLE( szFileName, *PI , myrank);
   
    READ_INT (szFileName, *iproc, myrank);
    READ_INT (szFileName, *jproc, myrank);
    
    *dx = *xlength / (double)(*imax);
    *dy = *ylength / (double)(*jmax);
    
    return 1;
}


/*Create the init_matrix for allocating the initial values and set the initial values*/
void init_uvp(
              double UI,
              double VI,
              double PI,
              int il,
              int ir,
              int jb,
              int jt,
              double **U,
              double **V,
              double **P
              ){
    
    init_matrix(U, il-2, ir+1, jb-1, jt+1, UI);
    init_matrix(V, il-1, ir+1, jb-2, jt+1, VI);
    init_matrix(P, il-1, ir+1, jb-1, jt+1, PI);
    
}



