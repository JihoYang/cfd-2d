/////////////////////////////////////////////////////////////////////////
//                                                                     //
//  Lattice Boltzmann Method (D3Q19) CFD Solver for Cavity Problem     //
//                                                                     //
//  Developer : Jiho Yang (MEng)                                       //
//              M.Sc. student, Computational Science & Engineering     //
//              Technische Universitat Munchen                         //
//                                                                     //
//  Final update date: 17/07/2016                                      //
//                                                                     //
/////////////////////////////////////////////////////////////////////////

#ifndef _MAIN_C_
#define _MAIN_C_

#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "LBDefinitions.h"
#include "helper.h"

int main(int argc, char *argv[]){

    //CPU time variables
    clock_t start_t, end_t, total_t;

    //Output directory variables
    char problem[50];
    char filename[50];
    char directory[50];
    char old_output_filename[128];
    struct dirent *old_outputfile;
    DIR  *output_directory;

    //Data fields
    double *collideField = NULL;
    double *streamField = NULL;
    int *flagField = NULL;
    int xlength;
    int numCells;

    //Initial conditions
    double velocityWall[3] = {0};		//Note that the velocityWall values (U_x, U_y, U_z) are initialised here 
    double velocityWall_X;			    //velocityWall in X direction 

    //Flow conditions
    double Re;
    double tau;
    double M;

    //Time variables
    int timesteps;
    int timestepsPerPlotting;
    int t;

    //Start CPU time measurement
    start_t = clock();

    //Read name of the problem from the command line (Please make sure that the name of the problem is the same as the .dat file)
    if(argc > 1){
 
        strcpy(problem, argv[1]);

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

    //Delete existing files in output folder
    while((old_outputfile = readdir(output_directory))){

        sprintf(old_output_filename, "%s/%s", problem, old_outputfile -> d_name);
        remove(old_output_filename);

    }
    
    printf("\n");

    //Read parameters
    readParameters(filename, &xlength, &velocityWall_X, &Re, &timesteps, &timestepsPerPlotting, argc, argv);

    printf("\n");
    printf("Simulation Start\n");
    printf("\n");

    //Set the first component of the array velocityWall to be the lid velocity in X direction 
    velocityWall[0] = velocityWall_X;	    

    //Compute Tau
    tau = (velocityWall_X * xlength)/(Re * pow(C_S, 2)) + 0.5;
    
    //Compute Mach number
    M = velocityWall_X / C_S;

    printf("\n");
    printf("Tau = %f    Mach = %f   Re = %f\n", tau, M, Re);
    printf("\n");
	
	//Allocate memory
	numCells = pow(xlength + 2, D);
	collideField = malloc(Q * numCells * sizeof(double));
	streamField = malloc(Q * numCells * sizeof(double)); 
	flagField = malloc(numCells * sizeof(int));
	
	//Set initial values of the fields 
	initialiseFields(collideField, streamField, flagField, xlength);
	
	for(t = 0; t < timesteps; t++){
	
		//Streaming step
		double *swap = NULL;
		
		doStreaming(collideField, streamField, flagField, xlength);
		swap = collideField;
		collideField = streamField;
		streamField = swap;

		//Collide step 
		doCollision(collideField, flagField, &tau, xlength);
		
		//Set boundary values 
		treatBoundary(collideField, flagField, velocityWall, xlength);

        //Export output
		if(t % timestepsPerPlotting == 0){

		    writeVtkOutput(collideField, flagField, directory, t, xlength);

		}

        printf("\n");
        printf("Timestep = %d finished\n", t);
		
	}

    // Close the output folder
    closedir(output_directory);
	
	//Free allocated memory
	free(collideField);
	free(streamField);
	free(flagField);

    //End CPU time measurement
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

	return 0;

}

#endif

