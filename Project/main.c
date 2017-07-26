//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  CFD Lab Summer Semester 2016                                            //
//                                                                          //
//  Project : Fluid-Structure Interaction for solving FIV problem           //
//                                                                          //
//  Group 8:   Mohammed Asif Chand                                          //
//             Jiho Yang                                                    //
//                                                                          //
//  Final update date: 14/07/2016                                           //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

//Convention for body index

//////////////////////////////////////////////////////////
//                                                      //
//   body 3        body 7                               //
//                                      body 11         //
//   body 2        body 6                               //
//                            body 8    body 10         //
//   body 1        body 5                               //
//                                      body 9          //
//   body 0        body 4                               //
//                                                      //
//////////////////////////////////////////////////////////

#include "helper.h"
#include "visual.h"
#include "uvp.h"
#include "boundary_val.h"
#include "boundary.h"
#include "init.h"
#include <stdio.h>
#include <string.h>
#include "sor.h"
#include "boundary.h"
#include "FSI.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>

int main(int argn, char** args){

    //Time variables
    clock_t start_t, end_t, total_t;

    //Fluid variables
    double **U, **V, **P, **F, **G, **RS;
    int **Flag, **pic_original = NULL;

    char problem[50];
    char filename[50];
    char directory[50];
    DIR  *output_directory;

    double Re, UI, VI, PI, GX, GY, t_end, xlength, ylength, dt, dx, dy, alpha, omg, omg_f, tau, eps, dt_value, dP;
    double res = 0, t = 0, n = 0;
    int imax, jmax, itermax, it;
    int wl, wr, wt, wb;

    //FSI variables
    int body_count, num_body;
    int *body_center_x, *body_center_y;
    int *disp_pixel, *disp_pixel_old, *net_disp_pixel;
    int it_f;
    double m, c, k;
    double *y, *v;
    double *force, *force_old;
    double res_f, *res_f_body;
    double **y_matrix;

    //Cut cell variables
    double *dy_upp, *dy_low;
    int cut_cell_flag;

    //Start time measurement
    start_t = clock();

    //Read name of the problem from the command line (Please make sure that the name of the problem is the same as the .dat file)
    if(argn > 1){

        strcpy(problem, args[1]);
    
    }

    else{

        printf("Please provide the name of the problem!\n e.g. Run ./sim <problem> for a problem based on problem.dat\n");

        return 1;

    }

    //Specify the input to be read based on the command line argument
    strcpy(filename, problem);
    strcat(filename, ".dat");

    //Create a problem specific directory for saving the output
    strcpy(directory, problem);
    strcat(directory, "/");
    strcat(directory, problem);
    mkdir(problem, 0777);
    output_directory = opendir(problem);

    //Read the program configuration file using read_parameters()
    read_parameters(filename, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &omg_f, &tau, &itermax, &eps, &dt_value,
                    &wl, &wr, &wt, &wb, &dP, &m, &c, &k, &num_body, &cut_cell_flag);

    printf("\n");
    printf("Simulation Start\n");
    printf("\n");

    //Set up the matrices (arrays) needed using the matrix() command
    U = matrix(0, imax  , 0, jmax+1);
    V = matrix(0, imax+1, 0, jmax  );
    P = matrix(0, imax+1, 0, jmax+1);
    F = matrix(0, imax  , 0, jmax+1);
    G = matrix(0, imax+1, 0, jmax  );
    RS= matrix(0, imax+1, 0, jmax+1);
    Flag = imatrix(0,imax+1,0,jmax+1);
    pic_original = imatrix(0, imax+1, 0, jmax+1);

    //Set up arrays for FSI variables
    force = malloc(num_body * sizeof(double));

    body_center_x = malloc(num_body * sizeof(int));
    body_center_y = malloc(num_body * sizeof(int));

    y = calloc(num_body, sizeof(double));
    v = calloc(num_body, sizeof(double));

    disp_pixel     = malloc(num_body * sizeof(int));
    disp_pixel_old = malloc(num_body * sizeof(int));
    net_disp_pixel = malloc(num_body * sizeof(int));

    res_f_body = malloc(num_body * sizeof(double));

    dy_upp = malloc(num_body*sizeof(double));
    dy_low = malloc(num_body*sizeof(double));
    
    y_matrix = matrix(0, t_end/dt, 0, num_body);
    

    //Set the initial conditions by using init_uvp and init_flag from init.c
    init_uvp(UI, VI, PI, imax, jmax, U, V, P);
    init_flag(problem, imax, jmax, Flag, pic_original);


    //Initialisation for displacement and velocity (for zero initial conditions, these lines are not necessary - calloc is used)
    for (body_count = 0; body_count <= num_body - 1; body_count++){

        y[body_count] = 0;
        v[body_count] = 0;
        dy_upp[body_count] = dy;
        dy_low[body_count] = dy;

    }

    int t_count = 0;

    //Iterate over time
    while(t <= t_end){

        it_f = 0;
        res_f = 1e6;
        force_old = calloc(num_body, sizeof(double));

        //Sub-iteration for force computation
        while (it_f < itermax && res_f > eps){

            //Compute the maximal time step size
            calculate_dt(Re, tau, &dt, dx, dy, dy_upp, dy_low, imax, jmax, U, V, num_body, cut_cell_flag);

            //Set the boundary conditions for each time step
            boundaryvalues(imax, jmax, U, V, wl, wr, wt, wb, Flag);

            //Set problem specific inflow boundary conditions
            spec_boundary_val(problem, imax, jmax, U, V, dP, Re, dx, dy, Flag);

            //Solve F(n) and G(n)
            calculate_fg(Re, GX, GY, alpha, dt, dx, dy, dy_upp, dy_low, imax, jmax, U, V, F, G, Flag, body_center_x, body_center_y, num_body, cut_cell_flag);

            //Solve the RHS of the Pressure Poisson Equation
            calculate_rs(dt, dx, dy, dy_upp, dy_low, imax, jmax, F, G, RS, Flag, body_center_x, body_center_y, num_body, cut_cell_flag);

            //Solve the whole Pressure Poisson Equation and compute P(n+1)
            it = 0;
            res = 1e6;

                while(it < itermax && res > eps){

                    sor(omg, dx, dy, dy_upp, dy_low, dP, imax, jmax, P, RS, &res, problem, Flag, body_center_x, body_center_y, num_body, cut_cell_flag);

                    it++;

                }

            //Update the velocities
            calculate_uv(dt, dx, dy, dy_upp, dy_low, imax, jmax, U, V, F, G, P, Flag, body_center_x, body_center_y, num_body, cut_cell_flag);

            //Calculate force
            calculate_force(P, force, dx, imax, jmax, Flag, num_body, body_center_x, body_center_y);

            //Calculate force residual
            for (body_count = 0; body_count <= num_body - 1; body_count++){

                //Under-relaxation for computing force
                force[body_count] = omg_f * force[body_count] + (1 - omg_f) * force_old[body_count];

                //Residual
                res_f_body[body_count] = fabs(force[body_count] - force_old[body_count]);
                force_old[body_count] = force[body_count];

                res_f += res_f_body[body_count];
                it_f++;

            }

            //Calculate average force residual
            res_f = res_f / num_body;

        }

        //Solve the equation of motion
        runge_kutta(y, v, m, c, k, force, dt, jmax, num_body);

        if( cut_cell_flag == CUT_CELL){
       
            //Cut the cell and convert the diplacement into pixels/grids
	        cut_cell(y, dx, dy, dy_upp, dy_low, imax, jmax, disp_pixel, disp_pixel_old, net_disp_pixel, num_body);

        }

        else{

	        //Convert displacement into number of pixels/grids
            disp2pixel(y, dx, dy, imax, jmax, disp_pixel, disp_pixel_old, net_disp_pixel, num_body);

        }

        //Update flag
        update_flag(imax, jmax, Flag, pic_original, net_disp_pixel, num_body, body_center_x, body_center_y);

        t = t + dt;
        n++;

        printf(" Timestep = %f | Residual (SOR) = %f | Iteration (SOR) = %d | Residual (force) = %f | Iteration (force) = %d\n", t, res, it, res_f, it_f);
        printf("\n");

        //Fill in y_matrix
        for (body_count = 1; body_count <= num_body; body_count++){

            y_matrix[t_count][body_count] = y[body_count-1];

        }
    
        y_matrix[t_count][0] = t;

        t_count++;

        //Export the solutions for visualisation
        if((int) n % (int) dt_value == 0){

            write_vtkFile(directory, n, xlength, ylength, imax, jmax, dx, dy, U, V, P);

        }


    }

    //Export displacement value over time
    export_disp(directory, dt, t, t_end, y_matrix);

    //Close the output folder
    closedir(output_directory);

    //Free memory allocation
    free_matrix(U , 0, imax  , 0, jmax+1);
    free_matrix(V , 0, imax+1, 0, jmax  );
    free_matrix(P , 0, imax+1, 0, jmax+1);
    free_matrix(F , 0, imax  , 0, jmax+1);
    free_matrix(G , 0, imax+1, 0, jmax  );
    free_matrix(RS, 0, imax+1, 0, jmax+1);
    free_matrix(y_matrix, 0, t_end/dt, 0, num_body);
    free_imatrix(Flag, 0, imax+1, 0, jmax+1);

    end_t = clock();
    total_t = (long double)(end_t - start_t) / CLOCKS_PER_SEC;

    printf("\n");
    printf("\n");
    printf("Total time taken by CPU: %lu\n", total_t);
    printf("\n");
    printf("Exiting of the program...\n");
    printf("\n");
    printf("Please find the output in the <%s> directory\n", problem);
    printf("\n");

    return -1;

}

