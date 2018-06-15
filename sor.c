#include <stdbool.h>
#include "sor.h"
#include "helper.h"
#include "boundary_val.h"

static short sorCounter = 0;

void computeSorSweep(double omg, double dx, double dy, int imax, int jmax, double **P, double **RS,
                     int **Flags, int iStart, int iEnd, int iStride, int jStart, int jEnd, int jStride);

void sor(double omg, double dx, double dy, int imax, int jmax, double **P, double **RS, int **Flags, double *res,
         int noFluidCells, double **U, double **V, short xFlowDirection, short yFlowDirection)
{
    double rloc;
    
    // This logic is to have SOR to alternate iteration directions
    short iSwitch = (short) ((sorCounter >> 1) & 1);
    short jSwitch = (short) (sorCounter & 1);
    int iStart = 1 + iSwitch*imax, iEnd = (!iSwitch)*imax + 1, iStride = 1-(2*iSwitch);
    int jStart = 1 + jSwitch*jmax, jEnd = (!jSwitch)*jmax + 1, jStride = 1-(2*jSwitch);
    /* SOR iteration */
    computeSorSweep(omg, dx, dy, imax, jmax, P, RS, Flags, iStart, iEnd, iStride, jStart, jEnd, jStride);
//    computeSorSweep(omg, dx, dy, imax, jmax, P, RS, Flags, 1, imax+1, 1, 1, jmax+1, 1);
//    computeSorSweep(omg, dx, dy, imax, jmax, P, RS, Flags, imax+1, 1, -1, 1, jmax+1, 1);
    
    // Update the counter
    ++sorCounter;
    sorCounter%=4;
    
    
    /* compute the residual */
    rloc = 0;
    for (int i = 1; i <= imax; i++)
    {
        for (int j = 1; j <= jmax; j++)
        {
            int cell = Flags[i][j];
            // proceed if fluid
            if (isFluid(cell))
            {
                int isRightFluid = isNeighbourFluid(cell, RIGHT) || (i == imax);
                int isLeftFluid = isNeighbourFluid(cell, LEFT) || (i == 1);
                int isTopFluid = isNeighbourFluid(cell, TOP) || (j == jmax);
                int isBottomFluid = isNeighbourFluid(cell, BOT) || (j == 1);
                double term =
                        (isRightFluid * (P[i + 1][j] - P[i][j]) - isLeftFluid * (P[i][j] - P[i - 1][j])) / (dx * dx)
                        + (isTopFluid * (P[i][j + 1] - P[i][j]) - isBottomFluid * (P[i][j] - P[i][j - 1])) / (dy * dy)
                        - RS[i][j];
                rloc += term *
                        term;
            }
        }
    }
    rloc = rloc / noFluidCells;
    rloc = sqrt(rloc);
    /* set residual */
    *res = rloc;
    
    /* set boundary values on the domain */
    for (int i = 1; i <= imax; i++)
    {
        // Here assuming PI=0
//        if (isInflow(Flags[i][0]))
//            P[i][0] = -0.5*U[i][0]*U[i][0];
//        else
            P[i][0] = P[i][1] * !(isOutflow(Flags[i][0]));
//        if (isInflow(Flags[i][jmax+1]))
//            P[i][jmax+1] = -0.5*U[jmax+1][0]*U[jmax+1][0];
//        else
            P[i][jmax + 1] = P[i][jmax] * !(isOutflow(Flags[i][jmax + 1]));
    }
    for (int j = 1; j <= jmax; j++)
    {
        // Here assuming PI=0
//        if (isInflow(Flags[0][j]))
//            P[0][j] = -0.5*U[0][j]*U[0][j];
//        else
            P[0][j] = P[1][j] * !(isOutflow(Flags[0][j]));
//        if (isInflow(Flags[imax+1][j]))
//            P[imax+1][j] = -0.5*U[imax+1][j]*U[imax+1][j];
//        else
            P[imax + 1][j] = P[imax][j] * !(isOutflow(Flags[imax + 1][j]));
    }
    
    /* set boundary values on obstacle interface */
    for (int i = 1; i <= imax; i++)
    {
        for (int j = 1; j <= jmax; j++)
        {
            setCenterBoundaryValues(imax, jmax, P, Flags, i, j, false);
        }
    }
}

void computeSorSweep(double omg, double dx, double dy, int imax, int jmax, double **P, double **RS,
                     int **Flags, int iStart, int iEnd, int iStride, int jStart, int jEnd, int jStride)
{
    for (int i = iStart; i != iEnd; i += iStride)
//    for (int j = jStart; j != jEnd; j += jStride)
    {
//        for (int i = iStart; i != iEnd; i += iStride)
        for (int j = jStart; j != jEnd; j += jStride)
        {
            int cell = Flags[i][j];
// proceed if fluid
            if (isFluid(cell))
            {
                int isRightFluid = isNeighbourFluid(cell, RIGHT) || (i == imax);
                int isLeftFluid = isNeighbourFluid(cell, LEFT) || (i == 1);
                int isTopFluid = isNeighbourFluid(cell, TOP) || (j == jmax);
                int isBottomFluid = isNeighbourFluid(cell, BOT) || (j == 1);
                double coeff = omg / (
                        (isRightFluid + isLeftFluid) / (dx * dx)
                        + (isTopFluid + isBottomFluid) / (dy * dy)
                );
                P[i][j] = (1.0 - omg) * P[i][j]
                          + coeff *
                            (
                                    (isRightFluid * P[i + 1][j] + isLeftFluid * P[i - 1][j]) / (dx * dx)
                                    + (isTopFluid * P[i][j + 1] + isBottomFluid * P[i][j - 1]) / (dy * dy)
                                    - RS[i][j]);
            }
        }
    }
}

