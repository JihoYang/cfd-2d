#ifndef __INIT_H_
#define __INIT_H_

#define CUT_CELL 1

#define C_F     16
#define C_B     0
#define B_O     8
#define B_W     4
#define B_S     2
#define B_N     1  
#define B_SO    10  
#define B_SW    6   
#define B_NO    9
#define B_NW    5

int read_parameters(
                    const char *szFileName,        // name of the file 
                    double *Re,                    // reynolds number   
                    double *UI,                    // velocity x-direction 
                    double *VI,                    // velocity y-direction 
                    double *PI,                    // pressure 
                    double *GX,                    // gravitation x-direction 
                    double *GY,                    // gravitation y-direction 
                    double *t_end,                 // end time 
                    double *xlength,               // length of the domain x-dir
                    double *ylength,               // length of the domain y-dir
                    double *dt,                    // time step
                    double *dx,                    // length of a cell x-dir
                    double *dy,                    // length of a cell y-dir
                    int  *imax,                    // number of cells x-direction
                    int  *jmax,                    // number of cells y-direction
                    double *alpha,                 // uppwind differencing factor
                    double *omg,                   // relaxation factor 
                    double *omg_f,                 // relaxation factor for force iteration
                    double *tau,                   // safety factor for time step
                    int  *itermax,                 // max. number of iterations  
                                                   // for pressure per time step 
                    double *eps,                   // accuracy bound for pressure
                    double *dt_value,              // time for output 
                    int *wl,                       // left boundary flag
                    int *wr,                       // right boundary flag
                    int* wt,                       // top boundary flag
                    int* wb ,                      // bottom boundary flag 
                    double *dP,                    // pressure difference
                    double *m,                     // mass
                    double *c,                     // damping coefficient
                    double *k,                     // stiffness
                    int *num_body,                 // number of bodies
		            int *cut_cell_flag             // binary flag to set and reset cut cell scheme
		           );


void init_uvp(
              double UI,
              double VI,
              double PI,
              int imax,
              int jmax,
              double **U,
              double **V,
              double **P
             );


void init_flag(
               char* problem,
               int imax,
               int jmax,
               int **Flag,
               int **pic_original
              );

#endif

