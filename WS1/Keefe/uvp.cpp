#include "uvp.h"	
#include <math.h>


void calculate_dt(
  double Re,
  double tau,
  double *dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V
){
	double u_max, v_max;
	for(int i=0; i<imax; i++){
		for(int j=0; j<jmax; j++){
			if abs(U[i][j]) > u_max
				u_max = abs(U[i][j]);
			if abs(V[i][j]) > v_max
				v_max = abs(U[i][j]);
		}
	}
	if(tau > 0){
		*dt = tau* min(Re/2*(1/(dx)^2 + 1/(dy)^2), dx/u_max, dy/v_max);
	}
}
