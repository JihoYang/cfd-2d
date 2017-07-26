#ifndef __UVP_H__
#define __UVP_H__

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
                 );


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
                 );


void calculate_dt(
                  double Re,
                  double tau,
                  double *dt,
                  double dx,
                  double dy,
	              double *dy_upp,
	              double *dy_low,
                  int    imax,
                  int    jmax,
                  double **U,
                  double **V,
                  int    num_body,
                  int    cut_cell_flag
                );


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
                );

#endif
