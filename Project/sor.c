#include "sor.h"
#include <math.h>
#include "init.h"
#include <string.h>
#include <stdio.h>

void sor(
         double omg,
         double dx,
         double dy,
         double *dy_upp,
         double *dy_low,
         double dP,
         int    imax,
         int    jmax,
         double **P,
         double **RS,
         double *res,
         char *problem,
         int **Flag,
         int *body_center_x,
         int *body_center_y,
         int cut_cell_flag,
         int num_body
        ){

    int i, j, a;
    int C_F_total = 0;
    double rloc;
    double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));
    
    int body_count = 0;
    int body_count_old = 0;

    //SOR iteration
    for(i = 1; i <= imax; i++) {

        for(j = 1; j<=jmax; j++) {

            //SOR iteration without cut-cell
            if(cut_cell_flag != CUT_CELL){

                if(Flag[i][j] > C_F) {

                    P[i][j] = (1.0 - omg) * P[i][j] + coeff *((P[i+1][j] + P[i-1][j]) / (dx*dx) 
                            + (P[i][j+1] + P[i][j-1]) / (dy*dy) - RS[i][j]);

                }

            }

            //SOR iteration for cut-cell
            else if (cut_cell_flag == CUT_CELL){

                if (i == body_center_x[body_count] && j == body_center_y[body_count]){
    
                    for (a = i - 1; a <= i + 1; a++){

                        P[a][j+2] = (1.0 - omg) * P[a][j+2] + omg/(2.0*(1.0/(dx*dx)+1.0/(dy_upp[body_count]*dy_upp[body_count]))) *((P[a+1][j+2] + P[a-1][j+2]) 
                                                / (dx*dx) + (P[a][j+3] + P[a][j+1]) / (dy_upp[body_count]*dy_upp[body_count]) - RS[a][j+2]);

                        P[a][j-2] = (1.0 - omg) * P[a][j-2] +  omg/(2.0*(1.0/(dx*dx)+1.0/(dy_low[body_count]*dy_low[body_count]))) *((P[a+1][j-2] + P[a-1][j-2]) 
                                                / (dx*dx) + (P[a][j] + P[a][j-3]) / (dy_low[body_count]*dy_low[body_count]) - RS[a][j-2]);

                    }

                    body_count ++;

                }
                
                if (body_count != 0){
                    
                    body_count_old = body_count - 1;
                    
                }

                else if (Flag[i][j] > C_F && j != body_center_y[body_count_old]+2 && j != body_center_y[body_count_old]-2){

                    P[i][j] = (1.0 - omg) * P[i][j] + coeff *((P[i+1][j] + P[i-1][j]) / (dx*dx) + (P[i][j+1] + P[i][j-1]) / (dy*dy) - RS[i][j]);

                }
    
            }

        }

    }

    //Compute the residual
    rloc = 0;
    body_count = 0;
    body_count_old = 0;

    for(i = 1; i <= imax; i++){

        for(j = 1; j <= jmax; j++){

            //Residual computation without cut-cell
            if(cut_cell_flag != CUT_CELL){

                if(Flag[i][j] > C_F){

                    rloc +=((P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx)
                         + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])
                         * ((P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);

                    C_F_total++;

                }

            }
   
            //Residual computation with cut-cell 
            else if(cut_cell_flag == CUT_CELL){

                if (i == body_center_x[body_count] && j == body_center_y[body_count]){

                    for (a = i - 1; a <= i + 1; a++){

                        rloc +=  (((P[a+1][j+2]-2.0*P[a][j+2]+P[a-1][j+2])/(dx*dx)
                                + ( P[a][j+3]-2.0*P[a][j+2]+P[a][j+1])/(dy_upp[body_count]*dy_upp[body_count]) - RS[a][j+2])
                                * ((P[a+1][j+2]-2.0*P[a][j+2]+P[a-1][j+2])/(dx*dx) + (P[a][j+3]-2.0*P[a][j+2]+P[a][j+1])/(dy_upp[body_count]*dy_upp[body_count]) - RS[a][j+2]))
                                + (((P[a+1][j-2]-2.0*P[a][j-2]+P[a-1][j-2])/(dx*dx)
                                + ( P[a][j-1]-2.0*P[a][j-2]+P[a][j-3])/(dy_low[body_count]*dy_low[body_count]) - RS[a][j-2])
                                * ((P[a+1][j-2]-2.0*P[a][j-2]+P[a-1][j-2])/(dx*dx) + (P[a][j-1]-2.0*P[a][j-2]+P[a][j-3])/(dy_low[body_count]*dy_low[body_count]) - RS[a][j-2]));

                    }

                    body_count ++;

                }
                
                if (body_count != 0){
                    
                    body_count_old = body_count - 1;
                    
                }

                else if(Flag[i][j] > C_F && j != body_center_y[body_count_old]+2 && j != body_center_y[body_count_old]-2){

                    rloc +=   ((P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx)
                            + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])
                            * ((P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);

                }

            }

            C_F_total++;

        }

    }

    rloc = rloc / C_F_total;
    rloc = sqrt(rloc);

    //Set residual
    *res = rloc;

    // Set values for horizontal boundaries (Neumann BC)
    for(i = 1; i <= imax; i++) {

        P[i][0] = P[i][1];
        P[i][jmax+1] = P[i][jmax];

    }

    //Neumann pressure BC
    if(dP != 0){

        for(j = 0; j <= jmax; j++) {

            P[0][j] = 2 * dP - P[1][j];
            P[imax+1][j] = -P[imax][j];

        }

    }

    //Dirichlet pressure BC
    else{

        for(j = 1; j <= jmax; j++) {

            P[0][j] = P[1][j];
            P[imax+1][j] = P[imax][j];

        }

    }

    //BC for inner obstacle cells
    for(i = 1; i <= imax; i++) {

        for(j = 1; j<=jmax; j++) {

            switch(Flag[i][j]) {

                case B_N:

                    P[i][j] = P[i][j+1];

                    break;

                case B_S:

                    P[i][j] = P[i][j-1];

                    break;

                case B_O:

                    P[i][j] = P[i+1][j];

                    break;

                case B_W:

                    P[i][j] = P[i-1][j];

                    break;

                case B_NW:

                    P[i][j] = 0.5 * (P[i-1][j] + P[i][j+1]);

                    break;

                case B_NO:

                    P[i][j] = 0.5 * (P[i+1][j] + P[i][j+1]);

                    break;

                case B_SW:

                    P[i][j] = 0.5 * (P[i-1][j] + P[i][j-1]);

                    break;

                case B_SO:

                    P[i][j] = 0.5 * (P[i+1][j] + P[i][j-1]);

                    break;

            }

        }

    }

}

