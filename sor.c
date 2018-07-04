#include <stdbool.h>
#include "sor.h"
#include "helper.h"
#include "boundary_val.h"
#include "logger.h"

static short k = 2;

void sor(double omg, double dx, double dy, int imax, int jmax, double **P, double **RS, int **Flags,
         BoundaryInfo *boundaryInfo, double *res, int noFluidCells)
{
    double rloc;
    
    /* SOR iteration */
    int istart, jstart, ibound, jbound, iflip, jflip, icounter, jcounter;
    
    icounter = k&1;
    jcounter = (k>>1)&1;
    
    istart = (imax-1) * icounter + 1;
    jstart = (jmax-1) * jcounter + 1;
    
    ibound = (imax + 1) * !icounter;
    jbound = (jmax + 1) * !jcounter;
    
    iflip = pow(-1, icounter);
    jflip = pow(-1, jcounter);
    for (int i = istart; i*iflip < ibound; i += iflip)
    {
        for(int j = jstart; j*jflip < jbound; j += jflip)
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
    ++k;
    k %= 4;
    
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
                rloc += (term * term);
            }
        }
    }
    rloc = rloc / noFluidCells;
    rloc = sqrt(rloc);
    
    /* set residual */
    *res = rloc;
    
    /* set boundary values on the domain */
    setPressureOuterBoundaryValues(imax, jmax, P, Flags, boundaryInfo);
    
    /* set boundary values on obstacle interface */
    for (int i = 1; i <= imax; i++)
    {
        for (int j = 1; j <= jmax; j++)
        {
            setCenterBoundaryValues(imax, jmax, P, Flags, i, j, false);
        }
    }
}


