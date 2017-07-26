#include "helper.h"
#include "init.h"

int read_parameters(
                    const char *szFileName,        //filename
                    double *Re,                    //Reynolds number
                    double *UI,                    //Initial x velocity
                    double *VI,                    //Initial y velocity
                    double *PI,                    //Initial pressure
                    double *GX,                    //Gravity in x direction
                    double *GY,                    //Gravity in y direction
                    double *t_end,                 //End time
                    double *xlength,               //Length of domain in x direction
                    double *ylength,               //Length of domain in y direction
                    double *dt,                    //Time step
                    double *dx,                    //Length of grid size in x direction
                    double *dy,                    //Length of grid size in y direction
                    int  *imax,                    //Number of grid cells in x direction
                    int  *jmax,                    //Number of grid cells in y direction
                    double *alpha,                 //Constant for donor cell scheme
                    double *omg,                   //Relaxation factor for SOR
                    double *tau,                   //Safety factor for adaptive global time stepping scheme
                    int  *itermax,                 //Maximum number of iterations for SOR
                    double *eps,                   //Machine precision
                    double *dt_value,              //Initial time step
                    int *wl,                       //Flag for left wall
                    int *wr,                       //Flag for right wall
                    int *wt,                       //Flag for upper wall
                    int *wb,                        //Flag for bottom wall
                    double *dP                         //Pressure difference
                    ){

   READ_DOUBLE( szFileName, *xlength );
   READ_DOUBLE( szFileName, *ylength );

   READ_DOUBLE( szFileName, *Re    );
   READ_DOUBLE( szFileName, *t_end );
   READ_DOUBLE( szFileName, *dt    );

   READ_INT   ( szFileName, *imax );
   READ_INT   ( szFileName, *jmax );

   READ_DOUBLE( szFileName, *omg   );
   READ_DOUBLE( szFileName, *eps   );
   READ_DOUBLE( szFileName, *tau   );
   READ_DOUBLE( szFileName, *alpha );

   READ_INT   ( szFileName, *itermax );
   READ_DOUBLE( szFileName, *dt_value );

   READ_DOUBLE( szFileName, *UI );
   READ_DOUBLE( szFileName, *VI );
   READ_DOUBLE( szFileName, *GX );
   READ_DOUBLE( szFileName, *GY );
   READ_DOUBLE( szFileName, *PI );
   READ_DOUBLE( szFileName, *dP);
   READ_INT(szFileName, *wl);
   READ_INT(szFileName, *wr);
   READ_INT(szFileName, *wt);
   READ_INT(szFileName, *wb);

   *dx = *xlength / (double)(*imax);
   *dy = *ylength / (double)(*jmax);

   return 1;
}


/*Create the init_matrix for allocating the initial values and set the initial values*/
void init_uvp(
              double UI,
              double VI,
              double PI,
              int imax,
              int jmax,
              double **U,
              double **V,
              double **P
              ){

    init_matrix(U, 0, imax ,  0, jmax+1, UI);
    init_matrix(V, 0, imax+1, 0, jmax  , VI);
    init_matrix(P, 0, imax+1, 0, jmax+1, PI);

}

/* flag to initialise the condition */

void init_flag(char* problem, int imax, int jmax, int **Flag){

    /* pgm's files are acquired  in the problem definition case in boundary.c */
    /* then the flag's are developed based on the corresponding pgm file */
    /* output of pgm is a 2 x2 matrix with pixel information */
    char image_filename[50];
    int **pic;
    int i,j;
    
    strcpy(image_filename, problem);
    strcat(image_filename, ".pgm");
    pic = read_pgm(image_filename);

    /* setting the values of the pic to zeros and 1 */
    for(i=0; i<= imax;i++){

        for(j=0; j<= jmax; j++){

            if(pic[i][j] ==  255){

                pic[i][j] = 1;

            }

            else{

                pic[i][j] = 0;

            }

        }

    }

  /* pic stores the information on the bits of the pixels */
  /* in simple terminology 1 implies that the fluid cell(pixel) is a fluid cell and 0  implies the corresponding pixel at i,j is boundary cell */
  /* within the boundary cells there are different kinds of boundary cells as defined in macros of init.h */
  /* instead of rewrtiting the flag to C_B and then taking the conditions for different types of cells , I have directly taken the values of Flag to conclude different
  types of boundary cells thereby reducing rewriting and as well as an extended use of if loop */

  /* Flag converts the cells into legal ( as per the worksheet ) fluid type and and boundary type based on the code mentioned below */

    for(i =1; i<=imax; i++){

        for(j=1; j<=jmax; j++){

            Flag[i][j] =  (((( pic[i][j] << 1  | pic [i][j+1]) << 1 | pic[i][j-1]) << 1 | pic[i][j-1]) << 1 |  pic[i][j+1]) ;

                /* flag values of all the forbidden cells =  14 (01110), 7(00111), 13 (01101), 11(01011) */
                if ( (Flag[i][j] == 14) || (Flag[i][j] ==  7) || (Flag[i][j] == 13) || (Flag[i][j] == 11)){

                printf("Forbidden cell is found! Stop the execution\n");
                exit(0);

                }

        }

    }

}
