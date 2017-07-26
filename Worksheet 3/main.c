#include "helper.h"
#include "visual.h"
#include "uvp.h"
#include "boundary_val.h"
#include "boundary.h"
#include "init.h"
#include <stdio.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "sor.h"
#include "boundary.h"

int main(int argn, char** args){

    clock_t start_t, end_t, total_t;    

    /*Define the variables required*/
    double **U, **V, **P, **F, **G, **RS;
    int **Flag;

    char problem[50];
    char filename[50];
    char directory[50];
    DIR  *output_directory;

    double Re, UI, VI, PI, GX, GY, t_end, xlength, ylength, dt, dx, dy, alpha, omg, tau, eps, dt_value, dP;
    double res = 0, t = 0, n = 0;
    int imax, jmax, itermax, it;
    int k = 0;
    int wl, wr, wt, wb;

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


    /*Read the program configuration file using read_parameters()*/
    read_parameters(filename, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value,
                    &wl, &wr, &wt, &wb, &dP);    

    printf("\n");
    printf("Simulation Start\n");
    printf("\n");

    /*Set up the matrices (arrays) needed using the matrix() command*/
    U = matrix(0, imax  , 0, jmax+1);
    V = matrix(0, imax+1, 0, jmax  );
    P = matrix(0, imax+1, 0, jmax+1);
    F = matrix(0, imax  , 0, jmax+1);
    G = matrix(0, imax+1, 0, jmax  );
    RS= matrix(0, imax+1, 0, jmax+1);
    Flag = imatrix(0,imax+1,0,jmax+1); 
    
    /*Set the initial conditions by using init_uvp and init_flag from init.c*/
    init_uvp(UI, VI, PI, imax, jmax, U, V, P);
    init_flag(problem, imax, jmax, Flag);

    //Set further initialisation for flow over step
    if(strcmp(problem, "flow_over_step") == 0){

         init_matrix(U, 0, imax ,  0, jmax/2, 0);

     }
 

    //Iterate over time
    while(t <= t_end){
      
        /*Compute the maximal time step size*/
        calculate_dt(Re, tau, &dt, dx, dy, imax, jmax, U, V);
        
        /*Set the boundary conditions for each time step*/
        boundaryvalues(imax, jmax, U, V, wl, wr, wt, wb, Flag);
    
        //Set problem specific inflow boundary conditions
        spec_boundary_val(problem, imax, jmax, U, V, dP, Re, dx, dy);
        
        /*Solve F(n) and G(n)*/
        calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U, V, F, G, Flag);
        
        /*Solve the RHS of the Pressure Poisson Equation*/
        calculate_rs(dt, dx, dy, imax, jmax, F, G, RS);
        
        /*Solve the whole Pressure Poisson Equation and compute P(n+1)*/
        it = 0;
        res = 1e6;
        
            while(it < itermax && res > eps){
                
                sor(omg, dx, dy, dP, imax, jmax, P, RS, &res, problem, Flag);
                
                it++;
                
            }

        /*Update the velocities*/
        calculate_uv(dt, dx, dy, imax, jmax, U, V, F, G, P);

        //Calculate Force
        //calculate_force();

        //Solve the equation of motion
        //Use ceil() for displacement - pixel conversion
        
        //Update the flag field
        
        t = t + dt;
        n++;
        k++;

        printf(" Residual = %f     |      Timestep = %f\n", res, t);

        //Export the solutions for visualisation
        write_vtkFile(directory, n, xlength, ylength, imax, jmax, dx, dy, U, V, P); 

       }
    
    // Close the output folder
    closedir(output_directory);

    
    free_matrix(U , 0, imax  , 0, jmax+1);
    free_matrix(V , 0, imax+1, 0, jmax  );
    free_matrix(P , 0, imax+1, 0, jmax+1);
    free_matrix(F , 0, imax  , 0, jmax+1);
    free_matrix(G , 0, imax+1, 0, jmax  );
    free_matrix(RS, 0, imax+1, 0, jmax+1);
    free_imatrix(Flag, 0, imax+1, 0, jmax+1);

    end_t = clock();
    total_t = (long double)(end_t - start_t) / CLOCKS_PER_SEC;

    printf("\n");
    printf("\n");
    printf("Total time taken by CPU: %lu\n", total_t  );
    printf("\n");
    printf("Exiting of the program...\n");
    printf("\n");
    printf("Please find the output in the <%s> directory\n", problem);
    printf("\n");
    
    return -1;
    
}

