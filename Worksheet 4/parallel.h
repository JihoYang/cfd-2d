#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#define MASTER 0

#define     RIGHT   1
#define     LEFT    2
#define     UP      3
#define     DOWN    4

#define     VAR_U           1
#define     VAR_V           2
#define     VAR_P           3

void Program_Message(char *txt);
/* produces a stderr text output  */



void Programm_Sync(char *txt);
/* produces a stderr textoutput and synchronize all processes */



void Programm_Stop(char *txt);
/* all processes will produce a text output, be synchronized and finished */



void init_parallel (
                    int iproc,
                    int jproc,
                    int imax, 
                    int jmax,
                    int *myrank,
                    int *il,
                    int *ir,
                    int *jb,
                    int *jt,
                    int *rank_l,
                    int *rank_r,
                    int *rank_b,
                    int *rank_t,
                    int *omg_i,
                    int *omg_j,
                    int num_proc
                   );


int check_np (int *iproc, int *jproc, int num_proc, int myrank);

void correct_ijproc (int *iproc, int *jproc, int num_proc);

void send2buffer(
    double **matrix, 
    double *buffer, 
    int il, 
    int ir,
    int jb, 
    int jt
);

void recv2buffer(
    double **matrix,
    double *buffer,
    int il,
    int ir, 
    int jb, 
    int jt
);

void communicate(
    double **matrix,
    int il, 
    int ir, 
    int jb, 
    int jt,   
    int rank_l, 
    int rank_r, 
    int rank_b, 
    int rank_t,
    int direction,
    int variable,
    double *bufSend, double *bufRecv
);


void pressure_comm(
                    double **P,                                                                                                                                                  
                    int  il,
                    int ir,
                    int jb,
                    int jt,
                    int rank_l,
                    int rank_r,
                    int rank_b,
                    int rank_t,
                    double *bufSend,
                    double *bufRecv
                    //MPI_Status *status - not sure where I could use this
                    //int chunk  - not sure what this is meant to be doing?
);


void uv_comm(
             double **U,
             double **V,
             int il, 
             int ir,
             int jb, 
             int jt,
             int rank_l, 
             int rank_r,
             int rank_b, 
             int rank_t,
             double *bufSend,
             double *bufRecv 
);
      

 
 



