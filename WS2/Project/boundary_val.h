#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__


/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(int omg_i, int omg_j, int imax_local, int jmax_local, double **U, double **V);
void boundaryvalues_FG(int omg_i, int omg_j, int imax_local, int jmax_local, double **F, double **G, double **U, double **V);
void boundaryvalues_P(int omg_i, int omg_j, int imax_local, int jmax_local, double **P);

#endif
