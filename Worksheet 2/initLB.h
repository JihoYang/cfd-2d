#ifndef _INITLB_H_
#define _INITLB_H_
#include "helper.h"


//Reads the parameters for the lid driven cavity scenario from a config file
int readParameters(
                   const char *szFileName,		       //Reads the user-defined file name. Parameter name : "szFileName"
                   int *xlength,                       //Reads domain size. Parameter name: "xlength" 
                   double *velocityWall,               //Velocity of the lid. Parameter name: "characteristicvelocity" 
                   double *Re,				           //Reynolds number. Parameter name: "Re" 
                   int *timesteps,                     //Number of timesteps. Parameter name: "timesteps" 
                   int *timestepsPerPlotting,          //Timesteps between subsequent VTK plots. Parameter name: "vtkoutput" 
                   int argc,                           //Number of arguments. Should equal 2 (program + name of config file 
                   char *argv[]                        //argv[1] shall contain the path to the config file 
);

//Initialises the particle distribution functions and the flagfield
void initialiseFields(
                      double *collideField, 
                      double *streamField,
                      int *flagField, 
                      int xlength
);

#endif

