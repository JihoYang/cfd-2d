#include "sor.h"
#include <math.h>
#include <mpi.h>
#include "parallel.h"
#include "helper.h"

void sor(
         double omg,
         double dx,
         double dy,
         int    il,
         int    ir,
         int    jb,
         int    jt,
         int    imax,
         int    jmax,
         int    rank_l,
         int    rank_r,
         int    rank_b,
         int    rank_t,
         double **P,
         double **RS,
         double *res,
         double *bufSend,
         double *bufRecv
        ){

    int i,j;
    double rloc;
    double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));     

    // BC for lower wall
    if(jb == 1){

        for(i = il; i <= ir; i++){

        P[i][0] = P[i][1];

        }

    }

    // BC for upper wall 
    if(jt == jmax){

        for(i = il; i <= ir; i++){

        P[i][jmax+1] = P[i][jmax];

        }

    }

    // BC for left wall
    if(il == 1){
    
        for(j = jb; j <= jt; j++){
    
        P[0][j] = P[1][j]; 
    
        }

    }

    // BC for right wall
    if(ir == imax){
    
        for(j = jb; j <= jt; j++){

        P[imax+1][j] = P[imax][j];
    
        }

    }

    /* SOR iteration */
    for(i = il; i <= ir; i++){

        for(j = jb; j<=jt; j++){

            P[i][j] = (1.0-omg)*P[i][j]
                    + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);

        }

    }

    //Exchange pressure boundary strips
    pressure_comm(P, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, bufSend, bufRecv);

    /* compute the residual */
    rloc = 0;

    for(i = il; i <= ir; i++){

        for(j = jb; j <= jt; j++){

            rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
                    ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);

        }

    }

    //Set residual
    rloc = rloc/(imax*jmax);
    rloc = sqrt(rloc);
    
    //Find global residual
    MPI_Allreduce(&rloc, res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

}

