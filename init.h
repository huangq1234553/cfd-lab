#ifndef __INIT_H_
#define __INIT_H_

#include "boundary_val.h"

typedef enum RunningMode { COMPACT=0, EXTENDED=1 } RunningMode;

static char *const NULLSTRING = "NULLSTRING";

char isStringDefault(char* string);

/**
 * This operation initializes all the local variables reading a configuration
 * file. For every variable a macro like READ_INT() is called passing it the
 * input filename and the variable to be written to. This macro calls
 * an operation read_int() augmenting the parameter set with the name of the
 * variable to be read. The read_int() operation parses the input file, extracts
 * the value of the variable, sets the variable and finally prints some debug
 * information. This is possible as the macro adds the name of the variable to
 * be set. All the helper operations can be found within helper.h and helper.c.
 *
 * @param szFileName char pointer to the filename
 * @param Re         Reynolds number
 * @param UI         initial velocity in  x-direction - used by init_uvp()
 * @param VI         initial velocity y-direction - used by init_upv()
 * @param PI         initial pressure - used by init_upv()
 * @param GX         gravitation x-direction
 * @param GY         gravitation y-direction
 * @param t_end      end time (not discrete in time steps)
 * @param xlength    domain length x-direction
 * @param ylength    domain lenght y-direction
 * @param dt         time step length: dividing t_end by dt gives the number of
 *                   time steps to perform. Actually dt is determined by a
 *                   function, so manipulating this value within the 
 *                   configuration file should not affect the solution process
 *                   at all
 * @param dx         cell length x-direction
 * @param dy         cell length y-direction
 * @param imax       number of cells in x-direction
 * @param jmax       number of cells in Y-direction
 * @param alpha      uppwind-differencing-factor alpha
 * @param omg        relaxation factor omega
 * @param tau        safety parameter for time step calculation
 * @param itermax    max. number of pressure iterations
 * @param eps        tolerance limit for pressure calculation
 * @param dt_value   time steps for output (after how many time steps one should
 *                   write into the output file)
 * @param problem    the problem short string (no spaces please!)
 * @param geometry   /path/to/geometry.pgm file
 */
int read_parameters(const char *szFileName, double *Re, double *UI, double *VI, double *PI, double *GX, double *GY,
                    double *t_end, double *xlength, double *ylength, double *dt, double *dx, double *dy, int *imax,
                    int *jmax, double *alpha, double *omg, double *tau, int *itermax, int *itermaxPGM,
                    int *sorIterationsThreshold, double *eps, double *dt_value, char *problem, char *geometry,
                    char *geometryMask, BoundaryInfo *boundaryInfo, double *beta, double *TI, double *T_h, double *T_c,
                    double *Pr, double *x_origin, double *y_origin, double *minVelocity, double *maxVelocity,
                    double *percentPressure, double *percentVelocity, int *checkPressure, int *checkVelocity,
                    int *checkVortex, int *vortexSizeThreshold, double *vortexStrengthThreshold,
                    double *obstacleBudgetFraction, char *precice_config, char *participant_name, char *mesh_name,
                    char *read_data_name, char *write_data_name, int *adaptiveGeometryEnabled);

//void read_boundary_parameters_compact_mode(const char *szFileName, BoundaryInfo *boundaryInfo, double dx, double dy);
void read_boundary_parameters_extended_mode(const char *szFileName, BoundaryInfo *boundaryInfo, double dx, double dy,
                                            int imax,
                                            int jmax, char *geometryFileName);

void configureBoundary(BoundaryInfo *boundaryInfo, BoundarySide boundarySide, double dirichletU, double dirichletV,
                       double pressure, double temp, double qN, double k, double h);

/**
 * The arrays U,V and P are initialized to the constant values UI, VI and PI on
 * the whole domain.
 */
void init_uvpt(double UI, double VI, double PI, double TI, int imax, int jmax, double **U, double **V, double **P,
               double **T, int **Flags);

void init_flag(char *problem, char *geometry, char *mask, int imax, int jmax, int **Flag, int *fluidCellsCounter,
               int *couplingCellsCounter, RunningMode runningMode);


#endif

