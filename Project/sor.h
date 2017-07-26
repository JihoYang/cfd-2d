#ifndef __SOR_H_
#define __SOR_H_

void sor(
         double omg,
         double dx,
         double dy,
     	 double *dy_upp,
         double *dy_low,
         double dP,
         int    imax,
         int    jmax,
         double **P,
         double **RS,
         double *res,
         char   *problem,
         int    **Flag,
         int    *body_center_x,
         int    *body_center_y,
         int    cut_cell_flag,
         int    num_body
        );


#endif
