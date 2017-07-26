#include "FSI.h"
#include "init.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "helper.h"

void calculate_force(
                     double **P, 
                     double *force, 
                     double dx, 
                     int imax, 
                     int jmax, 
                     int **Flag, 
                     int num_body, 
                     int *body_center_x, 
                     int *body_center_y
                    ){

    int i, j;
    int body_count = 0;

    double *P_upp_sum;
    double *P_low_sum;

    double *P_upp_avg;
    double *P_low_avg;

    double *P_diff;

    //Allocate memory
    P_upp_sum = malloc(num_body * sizeof(double));
    P_low_sum = malloc(num_body * sizeof(double));
    P_upp_avg = malloc(num_body * sizeof(double));
    P_low_avg = malloc(num_body * sizeof(double));
    P_diff    = malloc(num_body * sizeof(double));

    //Calculate upper and lower pressure
    for (i = 1; i <= imax; i++){

        for (j = 1; j <= jmax; j++){

            if (Flag[i][j] == 0){

                body_center_x[body_count] = i;
                body_center_y[body_count] = j;

                P_upp_sum[body_count] = P[i+1][j+1] + P[i][j+1] + P[i-1][j+1];
                P_low_sum[body_count] = P[i+1][j-1] + P[i][j-1] + P[i-1][j-1];

                body_count ++;

            }

        }

    }

    //Calculate pressure difference and resulting force
    for (body_count = 0; body_count <= num_body - 1; body_count++){

        P_upp_avg[body_count] = P_upp_sum[body_count] / 3;
        P_low_avg[body_count] = P_low_sum[body_count] / 3;

        //Note that the sign convention for force is positive if up, and negative if down
        P_diff[body_count] = P_low_avg[body_count] - P_upp_avg[body_count];

        force[body_count] = P_diff[body_count] * 3 * dx;

    }

}


void runge_kutta(
                 double *y,
                 double *v,
                 double m,
                 double c,
                 double k,
                 double *force,
                 double dt,
                 int jmax,
                 int num_body
                ){

    //y = Location in y direction, v = velocity in y direction, a = acceleration in y direction
	double fm, cm, km;
    double y2, y3, y4, v1, v2, v3, v4, a1, a2, a3, a4;
    int body_count;

    for (body_count = 0; body_count <= num_body - 1; body_count++){

	    fm = force[body_count] / m;
	    cm = c / m;
        km = k / m;

        //Initial conditions
        //y1 = y[body_count];
        v1 = v[body_count];

        //Runge-Kutta intermediate steps
        a1 = fm - cm * v[body_count] - km * y[body_count];

        y2 = y[body_count] + 0.5 * v1 * dt;
        v2 = v[body_count] + 0.5 * a1 * dt;
        a2 = fm - cm * v2 -km * y2;

        y3 = y[body_count] + 0.5 * v2 * dt;
        v3 = v[body_count] + 0.5 * a2 * dt;
        a3 = fm - cm * v3 - km * y3;

        y4 = y[body_count] + v3 * dt;
        v4 = v[body_count] + a3 * dt;
        a4 = fm - cm * v4 - km * y4;

        //Update final solution (displacement and velocity)
        y[body_count] = y[body_count] + dt / 6 * (v1 + 2 * v2 + 2 * v3 + v4);
        v[body_count] = v[body_count] + dt / 6 * (a1 + 2 * a2 + 2 * a3 + a4);

    }

}


//Converts the displacement into number of pixels to be moved
//Note that we are always shifting from the "original pic" (from the very initial state)
//This way, we can find the net displacement by simply using the location of the body (which is the output from runge_kutta)
void disp2pixel(
                double *y,
                double dx,
                double dy,
                int imax,
                int jmax,
                int *disp_pixel,
                int *disp_pixel_old,
                int *net_disp_pixel,
                int num_body
               ){

    double num_grid;
    int body_count;

    for (body_count = 0; body_count <= num_body - 1; body_count++){

        num_grid = y[body_count] / dy;

        if (y[body_count] > 0 && y[body_count] > 0.5 * dy){

            disp_pixel[body_count] = ceil( num_grid );

        }

        else if (y[body_count] < 0 && fabs(y[body_count]) > 0.5 * dy){

            disp_pixel[body_count] = floor( num_grid );

        }

        net_disp_pixel[body_count] = disp_pixel[body_count] - disp_pixel_old[body_count];
        disp_pixel_old[body_count] = disp_pixel[body_count];

    }

}


void cut_cell( 
              double *y,
              double dx,
              double dy,
              double *dy_upp,
              double *dy_low,
              int imax,
              int jmax,
              int *disp_pixel,
              int *disp_pixel_old,
              int *net_disp_pixel,
              int num_body 
             ){

    double num_grid;
    int body_count;

        for (body_count = 0; body_count <= num_body - 1; body_count++){

            if (fabs(y[body_count]) < dy && y[body_count] > 0){

                dy_upp[body_count] = dy - y[body_count];
                dy_low[body_count] = dy + y[body_count];

            }

            else if (fabs(y[body_count]) < dy && y[body_count] < 0){

                dy_upp[body_count] = dy + fabs(y[body_count]);
                dy_low[body_count] = dy - fabs(y[body_count]);

            }

            else if (fabs(y[body_count]) > dy && y[body_count] > 0){

                dy_upp[body_count] = dy - fmod(y[body_count], dy);
                dy_low[body_count] = dy + fmod(y[body_count], dy);

            }

            else if (fabs(y[body_count]) > dy && y[body_count] < 0){

                dy_upp[body_count] = dy + fmod(fabs(y[body_count]), dy);
                dy_low[body_count] = dy - fmod(fabs(y[body_count]), dy);

            }

            //Convert displacement (location to be precise) to number of pixels when criterion is met
            num_grid = y[body_count] / dy;

            if (y[body_count] > 0 && fmod(y[body_count], dy) > 0.9*dy){

                disp_pixel[body_count] = ceil( num_grid );

            }

            else if (y[body_count] < 0 && fmod(fabs(y[body_count]), dy) > 0.9*dy){

                disp_pixel[body_count] = floor( num_grid );

            }

            //Get net displacement (pixel)
            net_disp_pixel[body_count] = disp_pixel[body_count] - disp_pixel_old[body_count];
            disp_pixel_old[body_count] = disp_pixel[body_count];

    }

}

//Updates flagfield
void update_flag(
                 int imax,
                 int jmax,
                 int **Flag,
                 int **pic_original,
                 int *net_disp_pixel,
                 int num_body,
                 int *body_center_x,
                 int *body_center_y
                ){
    
    //disp_pixel is the displacement of the object in pixel (number of pixels to be moved)
    int i, j, a, b;
    int **pic_new;
    int body_count = 0;
    
    pic_new = imatrix(0,imax+1,0,jmax+1);
    init_imatrix(pic_new, 0, imax+1 ,  0, jmax+1, 1);
    
    //Create pic_new
    for (i = 1; i <= imax; i++){
        
        for (j = 1; j <= jmax; j++){
            
            if (i == body_center_x[body_count] && j == body_center_y[body_count]){
                
                for (a = i - 1; a <= i + 1; a++){
                    
                    for (b = j - 1; b <= j + 1; b++){
                        
                        pic_new[a][b + net_disp_pixel[body_count]] = 0;
                        
                    }
                    
                }
                
                body_count ++;
                
            }
            
        }
        
    }
    
    //Check if structure is out of domain
    if (body_count < num_body){
        
        printf("Structure out of domain - please adjust material properties\n");
        printf("\n");
        printf("Aborting program...\n");
        printf("\n");
        
        exit(0);
        
    }
    
    //Update Flag
    for(i = 1; i <= imax; i++){
        
        for(j = 1; j <= jmax; j++){
            
            Flag[i][j] =  (((( pic_new[i][j] << 1  | pic_new[i+1][j]) << 1 | pic_new[i-1][j]) << 1 | pic_new[i][j-1]) << 1 |  pic_new[i][j+1]);
            
        }
        
    }
    
}

