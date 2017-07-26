#ifndef _FSI_H_
#define _FSI_H_

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
                    );

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
                );

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
               );


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
             );


void update_flag(
                 int imax,
                 int jmax,
                 int **Flag,
                 int **pic_original,
                 int *net_disp_pixel,
                 int num_body,
                 int *body_center_x,
                 int *body_center_y
                );




#endif
