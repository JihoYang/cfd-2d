#ifndef __VISUAL_H__
#define __VISUAL_H__

void export_disp(                                                                                                                                    
                  const char *szProblem,
                  double dt,
                  double t,
                  double t_end,
                  double **y_matrix
                );


void write_vtkFile(
                   const char *szProblem,
		           int    timeStepNumber,
		           double xlength,
                   double ylength,
                   int    imax,
                   int    jmax,
		           double dx,
		           double dy,
                   double **U,
                   double **V,
                   double **P
                  );


void write_vtkHeader(
                     FILE *fp, 
                     int imax, 
                     int jmax, 
                     double dx, 
                     double dy
                    );


void write_vtkPointCoordinates(
                               FILE *fp, 
                               int imax, 
                               int jmax, 
                               double dx, 
                               double dy
                              );


#endif

