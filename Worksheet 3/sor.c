#include "sor.h"
#include <math.h>
#include "init.h"
#include <string.h>


void sor(
         double omg,
         double dx,
         double dy,
         double dP,
         int    imax,
         int    jmax,
         double **P,
         double **RS,
         double *res,
         char* problem,
         int** Flag
        ){

  int i,j;
  double rloc;
  double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));
  int C_F_total; /* total fluid cells */

  C_F_total = 0;

    /* SOR iteration */
    for(i = 1; i <= imax; i++) {

        for(j = 1; j<=jmax; j++) {

            if(Flag[i][j] == C_F){

            P[i][j] = (1.0-omg)*P[i][j] + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);

            }

        }

    }

    /* compute the residual */
    rloc = 0;

    for(i = 1; i <= imax; i++) {

        for(j = 1; j <= jmax; j++) {

            if( Flag[i][j] == C_F){   /* define all the cell C_f when proble */

            rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
                    ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);

       C_F_total++; /* total fluid cells */

            }

        }

    }


  rloc = rloc/(C_F_total);
  rloc = sqrt(rloc);

  /* set residual */
  *res = rloc;



    /* set pressure boundary conditions at the walls*/
    for(i = 1; i <= imax; i++){
           
        P[i][0] = P[i][1];
        P[i][jmax+1] = P[i][jmax];
           
    }

   for(j = 1; j <= jmax; j++) {
     
       P[0][j] = P[1][j];
       P[imax+1][j] = P[imax][j];
     
   }
      
    //IMPORTANT 
    /* defining boundary conditions for p */ //This is Neumann pressure BC which needs to be set for dP == 0
    for(i =0; i<= imax; i++){
        
        for( j=0; j<=jmax; j++){
            
            if(Flag[i][j] == B_O)
                
                P[i][j] = P[i+1][j];
            
            else if( Flag[i][j] == B_SO)
                
                P[i][j] = 0.5*(P[i+1][j] + P[i][j-1]) ;
            
            else if( Flag[i][j] == B_S)
                
                P[i][j] = P[i][j-1];
            
            else if( Flag[i][j] == B_SW)
                
                P[i][j] = 0.5*( P[i][j-1] + P[i-1][j]);
            
            else if( Flag[i][j] == B_NO)
                
                P[i][j] = 0.5* (P[i][j+1] + P[i+1][j]);
            
            else if( Flag[i][j] == B_W)
                
                P[i][j] = P[i-1][j];
            
            else if( Flag[i][j] == B_NW)
                
                P[i][j] = 0.5*(P[i][j+1] + P[i-1][j]);
            
            else if( Flag[i][j] == B_N)
                
                P[i][j] = P[i][j+1];
            
        }
        
    }
   


}
