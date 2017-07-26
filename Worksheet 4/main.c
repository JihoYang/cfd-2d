 //////////////////////////////////////////////////////////////////////////////
 //                                                                          //
 //  CFD Lab Summer Semester 2016                                            //
 //                                                                          //  
 //  Worksheet 4 : Parallelisation of 2D Finite Difference CFD solver        //
 //                                                                          // 
 //  Group 8:   Mohammed Asif Chand                                          //
 //             Jiho Yang                                                    //
 //                                                                          //
 //  Final update date: 15/06/2016                                           //
 //                                                                          //
 ////////////////////////////////////////////////////////////////////////////// 

//Please note that some modifications were made in helper.c (read_int and read_double) such that it includes myrank 
//This is purely for printing texts 

#include <stdio.h>
#include "helper.h"
#include "init.h"
#include "boundary_val.h"
#include "sor.h"
#include "uvp.h"
#include "visual.h"
#include "parallel.h"   

int main(int argn, char **args){
    
    //Define the variables required
    clock_t start_t, end_t, total_t; 
    double **U, **V, **P, **F, **G, **RS; 
    double *bufSend, *bufRecv;
    const char *szFileName = "cavity100.dat";
    double Re, UI, VI, PI, GX, GY, t_end, xlength, ylength, dt, dx, dy, alpha, omg, tau, eps, dt_value;
    double min_dt;
    double res = 0;
    double t = 0;
    double n = 0;
    int imax, jmax, itermax; 
    int it;

    //Variables required for parallelisation
    int iproc, jproc, myrank, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, omg_i, omg_j, num_proc;

    start_t = clock();

    //Initialise MPI environment
    MPI_Init(&argn, &args);

    //Get number of MPI tasks (np value) 
    //np is not number of processors being used, but number of processes which our program (which used to be one process if done sequentially) will be divided into
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc); 
    
    //Get number of rank (taskID of each process)
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 

    //Read the program configuration file using read_parameters()
    read_parameters(szFileName, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value, 
                    &iproc, &jproc, myrank);

    //Check if number of processes is appropriate and, if necessary, make changes accordingly
    check_np (&iproc, &jproc, num_proc, myrank);
    
    MPI_Barrier(MPI_COMM_WORLD);

    if (myrank == 0){

    printf("\n");
    printf("Decomposing domain\n");
    printf("\n");

    }

    //Domain decomposition and define neighbours
    init_parallel(iproc, jproc, imax, jmax, &myrank, &il, &ir, &jb, &jt, &rank_l, &rank_r, &rank_b, &rank_t, &omg_i, &omg_j, num_proc);
    
    //MPI_Barrier() is already included in init_parallel, but when having small number np it still somehow gets interrupted. Additional MPI_Barrier() is written here just to be sure... (although not quite sure why the MPI won't wait until all the other processes stop..!)
    MPI_Barrier(MPI_COMM_WORLD);

    if (myrank == 0){

    printf("\n");
    printf("Simulation Start\n");
    printf("\n");

    }
    
    //Set up the matrices (arrays) needed using the matrix() command
    U = matrix(il-2, ir+1, jb-1, jt+1);
    V = matrix(il-1, ir+1, jb-2, jt+1);
    P = matrix(il-1, ir+1, jb-1, jt+1);
    F = matrix(il-2, ir+1, jb-1, jt+1);
    G = matrix(il-1, ir+1, jb-2, jt+1);  
    RS= matrix(il, ir, jb, jt);

    
    //Set the initial conditions by using init_uvp from init.c
    init_uvp(UI, VI, PI, il, ir, jb, jt, U, V, P);

    //Create buffer space. This may be not necessary since (I believe) MPI manages the buffer memory space itself, but this makes the code more intuitive as of now..
    bufSend = malloc(max(ir - il + 2, jt - jb + 2) * sizeof(double));
    bufRecv = malloc(max(ir - il + 2, jt - jb + 2) * sizeof(double));

    while(t <= t_end){
        
        //Compute the maximal time step size
        calculate_dt(Re, tau, &dt, dx, dy, il, ir, jb, jt, U, V);

        //Find the minimum dt
        MPI_Allreduce(&dt, &min_dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        dt = min_dt;
        
        //Set the boundary conditions for each time step
        boundaryvalues(il, ir, jb, jt, imax, jmax, U, V);
        
        //Solve F(n) and G(n)
        calculate_fg(Re, GX, GY, alpha, dt, dx, dy, il, ir, jb, jt, imax, jmax, U, V, F, G);
        
        //Solve the RHS of the Pressure Poisson Equation
        calculate_rs(dt, dx, dy, il, ir, jb, jt, imax, jmax, F, G, RS);
    
        //Solve the whole Pressure Poisson Equation and compute P(n+1)
        it = 0;
        res = 1e6;
        
            while(it < itermax && res > eps){
                
                sor(omg, dx, dy, il, ir, jb, jt, imax, jmax, rank_l, rank_r, rank_b, rank_t, P, RS, &res, bufSend, bufRecv);
                
                it++;
                
            }

        //Update the velocities
        calculate_uv(dt, dx, dy, il, ir, jb, jt, imax, jmax, U, V, F, G, P);

        //Exchange velocities at boundary strips
        uv_comm(U, V, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, bufSend, bufRecv);
        
        t = t + dt;
        n++;

        if (myrank == 0){

        printf(" Residual = %f     |      Timestep = %f     |       Number of SOR Iteration = %d\n", res, t, it);

        }

        //Export the solutions for visualisation 
        if((int) n % (int) dt_value == 0){
            
        write_vtkFile("cavity", myrank, n, xlength, ylength, il, ir, jb, jt, imax, jmax, dx, dy, U, V, P);
            
        }
        
    }
    
    free_matrix(U, il-2, ir+1, jb-1, jt+1);
    free_matrix(V, il-1, ir+1, jb-2, jt+1);
    free_matrix(P, il-1, ir+1, jb-1, jt+1);
    free_matrix(F, il-2, ir+1, jb-1, jt+1);
    free_matrix(G, il-1, ir+1, jb-2, jt+1);
    free_matrix(RS, il, ir, jb, jt);
    free(bufSend);
    free(bufRecv);

    //Syncrhonise all processes
    MPI_Barrier(MPI_COMM_WORLD);    

    //End the MPI session
    MPI_Finalize(); 

    end_t = clock();
    total_t = (long double)(end_t - start_t) / CLOCKS_PER_SEC;

    if (myrank == 0){

        printf("\n");
        printf("\n");
        printf("Total time taken by CPU: %lu\n", total_t  );
        printf("\n");
        printf("Exiting of the program...\n");
        printf("\n");

    }
   
    return 0;
    
}

