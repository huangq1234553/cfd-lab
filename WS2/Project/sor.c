#include "sor.h"
#include "helper.h"
#include "boundary_val.h"

void sor(double omg, double dx, double dy, int imax, int jmax, double **P, double **RS, int **Flags, double *res, int noFluidCells)
{
    double rloc;
    double coeff = omg / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy)));
    
    /* SOR iteration */
    for (int i = 1; i <= imax; i++)
    {
        for (int j = 1; j <= jmax; j++)
        {
            int cell = Flags[i][j];
            // proceed if fluid
            if (isFluid(cell))
            {
                P[i][j] = (1.0 - omg) * P[i][j]
                          + coeff *
                            ((P[i + 1][j] + P[i - 1][j]) / (dx * dx) + (P[i][j + 1] + P[i][j - 1]) / (dy * dy) -
                             RS[i][j]);
            }
        }
    }

    
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
                rloc += ((P[i + 1][j] - 2.0 * P[i][j] + P[i - 1][j]) / (dx * dx) +
                         (P[i][j + 1] - 2.0 * P[i][j] + P[i][j - 1]) / (dy * dy) - RS[i][j]) *
                        ((P[i + 1][j] - 2.0 * P[i][j] + P[i - 1][j]) / (dx * dx) +
                         (P[i][j + 1] - 2.0 * P[i][j] + P[i][j - 1]) / (dy * dy) - RS[i][j]);
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
        P[i][0] = P[i][1];
        P[i][jmax + 1] = P[i][jmax];
    }
    for (int j = 1; j <= jmax; j++)
    {
        P[0][j] = P[1][j];
        P[imax + 1][j] = P[imax][j];
    }
    
    /* set boundary values on obstacle interface */
    for (int i = 1; i <= imax; i++)
    {
        for (int j = 1; j <= jmax; j++)
        {
            setCenterBoundaryValues(imax, jmax, P, Flags, i, j);
        }
    }
}

