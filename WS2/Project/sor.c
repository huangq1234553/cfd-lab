#include "sor.h"
#include "helper.h"
#include <math.h>

void sor(double omg, double dx, double dy, int imax, int jmax, double **P, double **RS, int **Flags, double *res, int noFluidCells)
{
    int i, j;
    double rloc;
    double coeff = omg / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy)));
    int cell;
    
//    /* set boundary values on obstacle interface */
    for (i = 1; i <= imax; i++)
    {
        for (j = 1; j <= jmax; j++)
        {
            cell = Flags[i][j];
            // proceed if obstacle
            if (isObstacle(cell))
            {
                if (isCorner(cell))
                {
                    P[i][j] = (P[i + isNeighbourObstacle(cell, LEFT) - isNeighbourObstacle(cell, RIGHT)][j] +
                               P[i][j + isNeighbourObstacle(cell, BOT) - isNeighbourObstacle(cell, TOP)]) / 2;
                }
                else
                {
                    P[i][j] = (isNeighbourFluid(cell, TOP)) * P[i][j + 1];
                    P[i][j] += (isNeighbourFluid(cell, BOT)) * P[i][j - 1];
                    P[i][j] += (isNeighbourFluid(cell, RIGHT)) * P[i + 1][j];
                    P[i][j] += (isNeighbourFluid(cell, LEFT)) * P[i - 1][j];
                }
            }
        }
    }
    
    /* SOR iteration */
    for (i = 1; i <= imax; i++)
    {
        for (j = 1; j <= jmax; j++)
        {
            cell = Flags[i][j];
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
    for (i = 1; i <= imax; i++)
    {
        for (j = 1; j <= jmax; j++)
        {
            cell = Flags[i][j];
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
    for (i = 1; i <= imax; i++)
    {
        P[i][0] = P[i][1];
        P[i][jmax + 1] = P[i][jmax];
    }
    for (j = 1; j <= jmax; j++)
    {
        P[0][j] = P[1][j];
        P[imax + 1][j] = P[imax][j];
    }

}

