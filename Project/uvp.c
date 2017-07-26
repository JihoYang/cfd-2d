#include "uvp.h"
#include <math.h>
#include <stdlib.h>
#include "helper.h"
#include "init.h"

/*Calculate the maximum time step allowed*/
void calculate_dt(
                  double Re,
                  double tau,
                  double *dt,
                  double dx,
                  double dy,
                  double *dy_upp,
                  double *dy_low,
                  int imax,
                  int jmax,
                  double **U,
                  double **V,
                  int num_body,
                  int cut_cell_flag
                 ){

    int i, j;
    double u_max = U[1][1], v_max = V[1][1];
    double dt_min = 0;
    double dy_min;

    int body_count;
    double dy_previous = 2*dy;

    //Find minimum dy among new dy_upp and dy_low values computed from cut-cell
    if(cut_cell_flag == CUT_CELL){

        for(body_count = 0 ; body_count < num_body-1; body_count++){

            dy_min = fmin(dy_upp[body_count], dy_previous);
            dy_min = fmin(dy_min, dy_low[body_count]);

            if(body_count != num_body -1){

                dy_previous = dy_min;
                
            }

            else{

                dy = fmin(dy_min,dy);
                
            }

        }
    
    }

    //Find the maximum velocities
    for (i = 1; i <= imax; i++){

        for (j = 1; j <= jmax; j++){

            if (fabs(U[i][j]) > u_max)

                u_max = fabs(U[i][j]);

            if (fabs(V[i][j]) > v_max)

                v_max = fabs(V[i][j]);

        }

    }

    //Time step criteron
    dt_min = fmin(dx / u_max, dy / v_max);
    dt_min = fmin(0.5 * Re / (1 / (dx * dx) + 1 / (dy * dy)), dt_min);

    if (tau < 0){

	    *dt = *dt;

    }

    else{

        *dt = tau * dt_min;

    }

}

//Calculate the F(n) and G(n)
void calculate_fg(
                  double Re,
                  double GX,
                  double GY,
                  double alpha,
                  double dt,
                  double dx,
                  double dy,
                  double *dy_upp,
                  double *dy_low,
                  int    imax,
                  int    jmax,
                  double **U,
                  double **V,
                  double **F,
                  double **G,
                  int    **Flag,
                  int    *body_center_x,
                  int    *body_center_y,
                  int    num_body,
                  int    cut_cell_flag
                 ){

    int i, j, a;
    int body_count = 0;
    int body_count_old = 0;

    //Boundary conditions for F (left and right)
    for(j = 1; j <= jmax; j++){

        F[0][j] = U[0][j];
        F[imax][j] = U[imax][j];

    }

    for (i = 1; i <= imax - 1; i++) {

        for (j = 1; j <= jmax; j++) {

            //Boundary conditions for inner cell obstacles
            if(Flag[i][j] == B_W) {

                F[i-1][j] = U[i-1][j];

            }

            else if(Flag[i][j] == B_O) {

                F[i][j] = U[i][j];

            }

            else if(Flag[i][j] == B_NO){

                F[i][j] = U[i][j];

            }

            else if(Flag[i][j] == B_NW){

                F[i-1][j] = U[i-1][j];

            }

            else if(Flag[i][j] == B_SO){

                F[i][j] = U[i][j];

            }

            else if(Flag[i][j] == B_SW){

                F[i-1][j] = U[i-1][j];

            }
    
            //Compute F without cut-cell
            if(cut_cell_flag != CUT_CELL){

                if(Flag[i][j] > C_F){

                    F[i][j] = U[i][j] + dt * (
                                           + 1 / Re * ((U[i+1][j] - 2 * U[i][j] + U[i-1][j]) / (dx * dx) + (U[i][j+1] - 2 * U[i][j] + U[i][j-1]) / (dy * dy))
                                           - 1 / dx * (pow((U[i][j] + U[i+1][j]) / 2, 2) - pow((U[i-1][j] + U[i][j]) / 2, 2))
                                           - alpha / dx * (fabs(U[i][j] + U[i+1][j]) * (U[i][j] - U[i+1][j]) / 4 - fabs(U[i-1][j] + U[i][j]) * (U[i-1][j] - U[i][j]) / 4)
                                           - 1 / dy * ((V[i][j] + V[i+1][j]) * (U[i][j] + U[i][j+1]) / 4 - (V[i][j-1] + V[i+1][j-1]) * (U[i][j-1] + U[i][j]) / 4)
                                           - alpha / dy * (fabs(V[i][j] + V[i+1][j]) * (U[i][j] - U[i][j+1]) / 4 - fabs(V[i][j-1] + V[i+1][j-1]) * (U[i][j-1] - U[i][j]) / 4)
                                           + GX
                                             );

                }

            }

            //Compute F with cut-cell
            else if(cut_cell_flag == CUT_CELL){

                if (i == body_center_x[body_count] && j == body_center_y[body_count]){

                    for (a = i - 1; a <= i + 1; a++){

                        F[a][j+2] = U[a][j+2] + dt * (
                                                   + 1 / Re * ((U[a+1][j+2] - 2 * U[a][j+2] + U[a-1][j+2]) 
                                                   / (dx * dx) + (U[a][j+3] - 2 * U[a][j+2] + U[a][j+1]) / (dy_upp[body_count] * dy_upp[body_count]))
                                                   - 1 / dx * (pow((U[a][j+2] + U[a+1][j+2]) / 2, 2) - pow((U[a-1][j+2] + U[a][j+2]) / 2, 2))
                                                   - alpha / dx * (fabs(U[a][j+2] + U[a+1][j+2]) * (U[a][j+2] - U[a+1][j+2]) / 4 
                                                   - fabs(U[a-1][j+2] + U[a][j+2]) * (U[a-1][j+2] - U[a][j+2]) / 4)
                                                   - 1 / dy_upp[body_count] * ((V[a][j+2] + V[a+1][j+2]) * (U[a][j+2] + U[a][j+3]) / 4 
                                                   - (V[a][j+1] + V[a+1][j+1]) * (U[a][j+1] + U[a][j+2]) / 4)
                                                   - alpha / dy_upp[body_count] * (fabs(V[a][j+2] + V[a+1][j+2]) * (U[a][j+2] - U[a][j+3]) / 4 
                                                   - fabs(V[a][j+1] + V[a+1][j+1]) * (U[a][j+1] - U[a][j+2]) / 4)
                                                   + GX
                                                     );

                        F[a][j-2] = U[a][j-2] + dt * (
                                                   + 1 / Re * ((U[a+1][j-2] - 2 * U[a][j-2] + U[a-1][j-2]) 
                                                   / (dx * dx) + (U[a][j-1] - 2 * U[a][j-2] + U[a][j-3]) / (dy_low[body_count] * dy_low[body_count]))
                                                   - 1 / dx * (pow((U[a][j-2] + U[a+1][j-2]) / 2, 2) - pow((U[a-1][j-2] + U[a][j-2]) / 2, 2))
                                                   - alpha / dx * (fabs(U[a][j-2] + U[a+1][j-2]) * (U[a][j-2] - U[a+1][j-2]) / 4 
                                                   - fabs(U[a-1][j-2] + U[a][j-2]) * (U[a-1][j-2] - U[a][j-2]) / 4)
                                                   - 1 / dy_low[body_count] * ((V[a][j-2] + V[a+1][j-2]) * (U[a][j-2] + U[a][j-1]) / 4 
                                                   - (V[a][j-3] + V[a+1][j-3]) * (U[a][j-3] + U[a][j-2]) / 4)
                                                   - alpha / dy_low[body_count] * (fabs(V[a][j-2] + V[a+1][j-2]) * (U[a][j-2] - U[a][j-1]) / 4 
                                                   - fabs(V[a][j-3] + V[a+1][j-3]) * (U[a][j-3] - U[a][j-2]) / 4)
                                                   + GX
                                                     );

                    }

                    body_count ++;

                }
                
                if (body_count != 0){
                    
                    body_count_old = body_count - 1;
                    
                }

                else if(Flag[i][j] > C_F && j != body_center_y[body_count_old]+2 && j != body_center_y[body_count_old]-2){

                    F[i][j] = U[i][j] + dt * (
                                           + 1 / Re * ((U[i+1][j] - 2 * U[i][j] + U[i-1][j]) / (dx * dx) + (U[i][j+1] - 2 * U[i][j] + U[i][j-1]) / (dy * dy))
                                           - 1 / dx * (pow((U[i][j] + U[i+1][j]) / 2, 2) - pow((U[i-1][j] + U[i][j]) / 2, 2))
                                           - alpha / dx * (fabs(U[i][j] + U[i+1][j]) * (U[i][j] - U[i+1][j]) / 4 - fabs(U[i-1][j] + U[i][j]) * (U[i-1][j] - U[i][j]) / 4)
                                           - 1 / dy * ((V[i][j] + V[i+1][j]) * (U[i][j] + U[i][j+1]) / 4 - (V[i][j-1] + V[i+1][j-1]) * (U[i][j-1] + U[i][j]) / 4)
                                           - alpha / dy * (fabs(V[i][j] + V[i+1][j]) * (U[i][j] - U[i][j+1]) / 4 - fabs(V[i][j-1] + V[i+1][j-1]) * (U[i][j-1] - U[i][j]) / 4)
                                           + GX
                                           );

                }

            }

        }

    }

    //Boundary conditions for G (bottom and top)
    for (i = 1; i <= imax; i++){

        G[i][0] = V[i][0];
        G[i][jmax] = V[i][jmax];

    }

    //Initialise body counter to zero for G computation
    body_count = 0;
    body_count_old = 0;

    for (i = 1; i <= imax; i++) {

        for (j = 1; j <= jmax - 1; j++) {

            //Boundary conditions for inner cell obstacles
            if(Flag[i][j] == B_N) {

                G[i][j] = V[i][j];

            }

            else if(Flag[i][j] == B_S) {

                G[i][j-1] = V[i][j-1];

            }

            else if(Flag[i][j] == B_NO) {

                G[i][j] = V[i][j];

            }

            else if(Flag[i][j] == B_NW) {

                G[i][j] = V[i][j];

            }

            else if(Flag[i][j] == B_SO) {

                G[i][j-1] = V[i][j-1];

            }

            else if(Flag[i][j] == B_SW) {

                G[i][j-1] = V[i][j-1];

            }

            //Compute G without cut-cell
            if(cut_cell_flag != CUT_CELL){

                if(Flag[i][j] > C_F) {

                    G[i][j] = V[i][j] + dt * (
                                           + 1 / Re * ((V[i][j+1] - 2 * V[i][j] + V[i][j-1]) / (dy * dy) + (V[i+1][j] - 2 * V[i][j] + V[i-1][j]) / (dx * dx))
                                           - 1 / dy * (pow((V[i][j] + V[i][j+1]) / 2, 2) - pow((V[i][j-1] + V[i][j]) / 2, 2))
                                           - alpha / dy * (fabs(V[i][j] + V[i][j+1]) * (V[i][j] - V[i][j+1]) / 4 - fabs(V[i][j-1] + V[i][j]) * (V[i][j-1] - V[i][j]) / 4)
                                           - 1 / dx * ((U[i][j] + U[i][j+1]) * (V[i][j] + V[i+1][j]) / 4 - (U[i-1][j] + U[i-1][j+1]) * (V[i-1][j] + V[i][j]) / 4)
                                           - alpha / dx * (fabs(U[i][j] + U[i][j+1]) * (V[i][j] - V[i+1][j]) / 4 - fabs(U[i-1][j] + U[i-1][j+1]) * (V[i-1][j] - V[i][j]) / 4)
                                           + GY
                                          );

                }

            }

            //Compute G with cut-cell
            else if(cut_cell_flag == CUT_CELL){
            
                if (i == body_center_x[body_count] && j == body_center_y[body_count]){

                    for (a = i - 1; a <= i + 1; a++){

                        G[a][j+2] = V[a][j+2] + dt * (
                                                   + 1 / Re * ((V[a][j+3] - 2 * V[a][j+2] + V[a][j+1]) / (dy_upp[body_count] * dy_upp[body_count]) 
                                                   + (V[a+1][j+2] - 2 * V[a][j+2] + V[a-1][j+2]) / (dx * dx))
                                                   - 1 / dy_upp[body_count] * (pow((V[a][j+2] + V[a][j+3]) / 2, 2) - pow((V[a][j+1] + V[a][j+2]) / 2, 2))
                                                   - alpha / dy_upp[body_count] * (fabs(V[a][j+2] + V[a][j+3]) * (V[a][j+2] - V[a][j+3]) / 4 
                                                   - fabs(V[a][j+1] + V[a][j+2]) * (V[a][j+1] - V[a][j+2]) / 4)
                                                   - 1 / dx * ((U[a][j+2] + U[a][j+3]) * (V[a][j+2] + V[a+1][j+2]) / 4 
                                                   - (U[a-1][j+2] + U[a-1][j+3]) * (V[a-1][j+2] + V[a][j+2]) / 4)
                                                   - alpha / dx * (fabs(U[a][j+2] + U[a][j+3]) * (V[a][j+2] - V[a+1][j+2]) / 4 
                                                   - fabs(U[a-1][j+2] + U[a-1][j+3]) * (V[a-1][j+2] - V[a][j+2]) / 4)
                                                   + GY
                                                     );

                        G[a][j-2] = V[a][j-2] + dt * (
                                                   + 1 / Re * ((V[a][j-1] - 2 * V[a][j-2] + V[a][j-3]) / (dy_low[body_count] * dy_low[body_count]) 
                                                   + (V[a+1][j-2] - 2 * V[a][j-2] + V[a-1][j-2]) / (dx * dx))
                                                   - 1 / dy_low [body_count]* (pow((V[a][j-2] + V[a][j-1]) / 2, 2) - pow((V[a][j-3] + V[a][j-2]) / 2, 2))
                                                   - alpha / dy_low[body_count] * (fabs(V[a][j-2] + V[a][j-1]) * (V[a][j-2] - V[a][j-1]) / 4 
                                                   - fabs(V[a][j-3] + V[a][j-2]) * (V[a][j-3] - V[a][j-2]) / 4)
                                                   - 1 / dx * ((U[a][j-2] + U[a][j-1]) * (V[a][j-2] + V[a+1][j-2]) / 4 
                                                   - (U[a-1][j-2] + U[a-1][j-1]) * (V[a-1][j-2] + V[a][j-2]) / 4)
                                                   - alpha / dx * (fabs(U[a][j-2] + U[a][j-1]) * (V[a][j-2] - V[a+1][j-2]) / 4 
                                                   - fabs(U[a-1][j-2] + U[a-1][j-1]) * (V[a-1][j-2] - V[a][j-2]) / 4)
                                                   + GY
                                                     );

                    }
        
                    body_count ++;

                }
                
                if (body_count != 0){
                    
                    body_count_old = body_count - 1;
                    
                }

                else if(Flag[i][j] > C_F && j != body_center_y[body_count_old]+2 && j !=  body_center_y[body_count_old]-2){

                    G[i][j] = V[i][j] + dt * (
                                           + 1 / Re * ((V[i][j+1] - 2 * V[i][j] + V[i][j-1]) / (dy * dy) + (V[i+1][j] - 2 * V[i][j] + V[i-1][j]) / (dx * dx))
                                           - 1 / dy * (pow((V[i][j] + V[i][j+1]) / 2, 2) - pow((V[i][j-1] + V[i][j]) / 2, 2))
                                           - alpha / dy * (fabs(V[i][j] + V[i][j+1]) * (V[i][j] - V[i][j+1]) / 4 - fabs(V[i][j-1] + V[i][j]) * (V[i][j-1] - V[i][j]) / 4)
                                           - 1 / dx * ((U[i][j] + U[i][j+1]) * (V[i][j] + V[i+1][j]) / 4 - (U[i-1][j] + U[i-1][j+1]) * (V[i-1][j] + V[i][j]) / 4)
                                           - alpha / dx * (fabs(U[i][j] + U[i][j+1]) * (V[i][j] - V[i+1][j]) / 4 - fabs(U[i-1][j] + U[i-1][j+1]) * (V[i-1][j] - V[i][j]) / 4)
                                           + GY
                                             );

                }

            }

        }

    }

}

//Calculate the RHS of the Pressure Poisson Equation
void calculate_rs(
                  double dt,
                  double dx,
                  double dy,
                  double *dy_upp,
                  double *dy_low,
                  int    imax,
                  int    jmax,
                  double **F,
                  double **G,
                  double **RS,
                  int    **Flag,
                  int    *body_center_x,
                  int    *body_center_y,
                  int    num_body,
                  int    cut_cell_flag
                 ){

    int i, j, a;
    int body_count = 0;
    int body_count_old = 0;

    for (i = 1; i <= imax; i++){

        for (j = 1; j <= jmax; j++){

            //Calculate residual without cut-cell
            if(cut_cell_flag != CUT_CELL){

                RS[i][j] = 1 / dt * ((F[i][j] - F[i-1][j]) / dx + (G[i][j] - G[i][j-1]) / dy);

            }

            //Calculate residual with cut-cell
            else if(cut_cell_flag == CUT_CELL){

                if (i == body_center_x[body_count] && j == body_center_y[body_count]){

                    for (a = i - 1; a <= i + 1; a++){

                        RS[a][j+2] = 1 / dt * ((F[a][j+2] - F[a-1][j+2]) / dx + (G[a][j+2] - G[a][j+1]) / dy_upp[body_count]);
                        RS[a][j-2] = 1 / dt * ((F[a][j-2] - F[a-1][j-2]) / dx + (G[a][j-2] - G[a][j-3]) / dy_low[body_count]);

                    }

                    body_count ++;

                }

                if (body_count != 0){
                    
                    body_count_old = body_count - 1;
                    
                }

                else if(Flag[i][j] > C_F && j != body_center_y[body_count_old]+2 && j != body_center_y[body_count_old]-2){

                    RS[i][j] = 1 / dt * ((F[i][j] - F[i-1][j]) / dx + (G[i][j] - G[i][j-1]) / dy);

                }

            }

        }

    }

}

//Update the velocities
void calculate_uv(
                  double dt,
                  double dx,
                  double dy,
                  double *dy_upp,
                  double *dy_low,
                  int    imax,
                  int    jmax,
                  double **U,
                  double **V,
                  double **F,
                  double **G,
                  double **P,
                  int    **Flag,
                  int    *body_center_x,
                  int    *body_center_y,
                  int    num_body,
                  int    cut_cell_flag
                 ){

    int i, j, a;
    int body_count = 0;
    int body_count_old = 0;

    for (i = 1; i <= imax - 1; i++){

        for (j = 1; j <= jmax; j++) {

            if(Flag[i][j] > C_F) {

                U[i][j] = F[i][j] - dt * (P[i+1][j] - P[i][j]) / dx;

            }

            else if(Flag[i][j] < C_F){

                U[i][j] = 0;

            }

        }

    }

    for (i = 1; i <= imax; i++) {

        for (j = 1; j <= jmax - 1; j++) {

            //Update UV without cut-cell
            if(cut_cell_flag != CUT_CELL){

                if(Flag[i][j] > C_F) {

                    V[i][j] = G[i][j] - dt * (P[i][j+1] - P[i][j]) / dy;

                }

                else{

                    V[i][j] = 0;

                }

            }

            //Update UV with cut-cell
            else if(cut_cell_flag == CUT_CELL){

                if (i == body_center_x[body_count] && j == body_center_y[body_count]){

                    for (a = i - 1; a <= i + 1; a++){

                        V[a][j+2] = G[a][j+2] - dt * (P[a][j+3] - P[a][j+2]) / dy_upp[body_count];
                        V[a][j-2] = G[a][j-2] - dt * (P[a][j] - P[a][j-2]) / dy_low[body_count];

                    }
        
                    body_count ++;

                }
                
                if (body_count != 0){
                    
                    body_count_old = body_count - 1;
                    
                }

                else if(Flag[i][j] > C_F && j != body_center_y[body_count_old]+2 && j != body_center_y[body_count_old]-2){

                    V[i][j] = G[i][j] - dt * (P[i][j+1] - P[i][j]) / dy;

                }

                else if (Flag[i][j] < C_F){

                    V[i][j] = 0;

                }

            }

        }

    }

}

