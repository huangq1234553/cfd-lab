#include "uvp.h"	
#include <math.h>


/**
 * Determines the value of F and G according to the formula
 *
 * @f$ F_{i,j} := u_{i,j} + \delta t \left( \frac{1}{Re} \left( \left[
    \frac{\partial^2 u}{\partial x^2} \right]_{i,j} + \left[
    \frac{\partial^2 u}{\partial y^2} \right]_{i,j} \right) - \left[
    \frac{\partial (u^2)}{\partial x} \right]_{i,j} - \left[
    \frac{\partial (uv)}{\partial y} \right]_{i,j} + g_x \right) @f$
 *
 * @f$ i=1,\ldots,imax-1, \quad j=1,\ldots,jmax @f$
 *
 * @f$ G_{i,j} := v_{i,j} + \delta t \left( \frac{1}{Re} \left(
   \left[ \frac{\partial^2 v}{\partial x^2}\right]_{i,j} + \left[ \frac{\partial^2 v}{\partial
                   y^2} \right]_{i,j} \right) - \left[ \frac{\partial
                   (uv)}{\partial x} \right]_{i,j} - \left[
                 \frac{\partial (v^2)}{\partial y} \right]_{i,j} + g_y
               \right) @f$
 *
 * @f$ i=1,\ldots,imax, \quad j=1,\ldots,jmax-1 @f$
 *
 */

void calculate_fg(double Re, double GX, double GY, double alpha, double dt, double dx, double dy, int imax, int jmax, double **U, double **V, double **F, double **G) {

	// set boundary conditions for F - see discrete momentum equations - apply Neumann BC - first derivative of pressure must be "zero" - dp/dx = 0
	for (int j = 1; j <= jmax; j++) {
		F[0][j] = U[0][j];
		F[imax][j] = U[imax][j];
	}

	// set boundary conditions for G - see discrete momentum equations - apply Neumann BC - first derivative of pressure must be "zero" - dp/dy = 0
	for (int i = 1; i <= imax; i++) {
		G[i][0] = V[i][0];
		G[i][jmax] = V[i][jmax];
	}

	// calculate F in the domain
	for (int i = 1; i < imax; i++) {
		for (int j = 1; j <= jmax; j++) {
			F[i][j] = ...
				// velocity u
				U[i][j] ...
				// diffusive term
				+ dt * ((1 / Re * ((U[i + 1][j] - 2 * U[i]][j] + U[i - 1][j]) / dx*dx + (U[i][j + 1] - 2 * U[i]][j] + U[i][j - 1]) / dy*dy)) ...
				// convective term
				- ( 1/dx * (((U[i][j] + U[i+1][j]) / 2)^2 - ((U[i-1][j] + U[i][j])/2)^2) + alpha/dx * ((abs(U[i][j] + U[i+1][j])/2 * (U[i][j] - U[i + 1][j])/2) - (abs(U[i-1][j] + U[i][j]) / 2 * (U[i - 1][j] - U[i][j]) / 2))) ...
				// convective term cont.
				- (1 / dy * ((V[i][j]+V[i+1][j])/2 * (U[i][j]+U[i][j+1])/2 - (V[i][j-1]+V[i+1][j-1]) / 2 * (U[i][j-1]+U[i][j]) / 2) + alpha / dy * (abs(V[i][j]+V[i+1][j])/2*(U[i][j] - U[i][j+1])/2 - abs(V[i][j-1]+V[i+1][j-1])/2 * (U[i][j-1] - U[i][j])/2)) ...
				// volume force
				+ GX[i][j]);
		}
	}

	// calculate G in the domain
	for (int i = 1; i <= imax; i++) {
		for (int j = 1; j < jmax; j++) {
			G[i][j] = ...
				// velocity v
				V[i][j] ...
				// diffusive term
				+ dt * ((1 / Re * ((V[i + 1][j] - 2 * V[i]][j] + V[i - 1][j]) / dx*dx + (V[i][j + 1] - 2 * V[i]][j] + V[i][j - 1]) / dy*dy)) ...
				// convective term
				- (1 / dy * (((V[i][j] + V[i][j + 1]) / 2) ^ 2 - ((V[i][j - 1] + V[i][j]) / 2) ^ 2) + alpha / dy * ((abs(V[i][j] + V[i][j + 1]) / 2 * (V[i][j] - V[i][j + 1]) / 2) - (abs(V[i][j - 1] + V[i][j]) / 2 * (V[i][j - 1] - V[i][j]) / 2))) ...
				// convective term cont.
				- (1 / dx * ((U[i][j] + U[i][j + 1]) / 2 * (V[i][j] + V[i + 1][j]) / 2 - (U[i - 1][j] + U[i - 1][j + 1]) / 2 * (V[i - 1][j] + V[i][j]) / 2) + alpha / dy * (abs(U[i][j] + U[i][j + 1]) / 2 * (V[i][j] - V[i + 1][j]) / 2 - abs(U[i - 1][j] + U[i - 1][j + 1]) / 2 * (V[i - 1][j] - V[i][j]) / 2)) ...
				// volume force
				+ GY[i][j]);
		}
	}
 }

/**
 * This operation computes the right hand side of the pressure poisson equation.
 * The right hand side is computed according to the formula
 *
 * @f$ rs = \frac{1}{\delta t} \left( \frac{F^{(n)}_{i,j}-F^{(n)}_{i-1,j}}{\delta x} + \frac{G^{(n)}_{i,j}-G^{(n)}_{i,j-1}}{\delta y} \right)  @f$
 *
 */
void calculate_rs(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **F,
  double **G,
  double **RS
){
  for(int i=0; i<imax+1; i++){
    for(int j=0; j<jmax+1; j++){
      RS[i][j] = ((F[i][j]-F[i-1][j])/dx + (G[i][j]-G[i][j-1])/dy)/dt;
    }
  }
}

/**
 * Determines the maximal time step size. The time step size is restricted
 * accordin to the CFL theorem. So the final time step size formula is given
 * by
 *
 * @f$ {\delta t} := \tau \, \min\left( \frac{Re}{2}\left(\frac{1}{{\delta x}^2} + \frac{1}{{\delta y}^2}\right)^{-1},  \frac{{\delta x}}{|u_{max}|},\frac{{\delta y}}{|v_{max}|} \right) @f$
 *
 */

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
		*dt = tau*min(Re/2*(1/pow(dx,2) + 1/pow(dy,2)), dx/u_max, dy/v_max);
	}
}

/**
 * Calculates the new velocity values according to the formula
 *
 * @f$ u_{i,j}^{(n+1)}  =  F_{i,j}^{(n)} - \frac{\delta t}{\delta x} (p_{i+1,j}^{(n+1)} - p_{i,j}^{(n+1)}) @f$
 * @f$ v_{i,j}^{(n+1)}  =  G_{i,j}^{(n)} - \frac{\delta t}{\delta y} (p_{i,j+1}^{(n+1)} - p_{i,j}^{(n+1)}) @f$
 *
 * As always the index range is
 *
 * @f$ i=1,\ldots,imax-1, \quad j=1,\ldots,jmax @f$
 * @f$ i=1,\ldots,imax, \quad j=1,\ldots,jmax-1 @f$
 *
 * @image html calculate_uv.jpg
 */

void calculate_uv(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G,
  double **P
)
{
  // TODO: [Tom] Check if here we must start from 0 or not (same for i/j-max) - boundary should not be updated...
  for (int i=1; i<imax; ++i)
  {
    for (int j=1; j<jmax; ++j)
    {
      U[i][j] = F[i][j] - dt*(P[i+1][j] - P[i][j])/dx;
      V[i][j] = G[i][j] - dt*(P[i][j+1] - P[i][j])/dy;
    }
  }
}