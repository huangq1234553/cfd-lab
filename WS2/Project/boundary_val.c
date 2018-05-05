#include "boundary_val.h"
#include "helper.h"

void boundaryvalues(int imax, int jmax, double **U, double **V, int **Flags, BoundaryInfo boundaryInfo[4])
{
    // Setting boundary conditions on the outer boundary
    leftboundary(imax, jmax, U, V, boundaryInfo);
    rightboundary(imax, jmax, U, V, boundaryInfo);
    topboundary(imax, jmax, U, V, boundaryInfo);
    bottomboundary(imax, jmax, U, V, boundaryInfo);
    
    // Boundary values at geometries in the internal part of the domain
    for (int i = 1; i <= imax; ++i)
    {
        for (int j = 1; j <= jmax; ++j)
        {
            int cell = Flags[i][j];
            if (isObstacle(cell))
            {
                // Compute v
                if (!skipV(cell))
                {
                    if (isNeighbourFluid(cell, TOP))
                    {
                        V[i][j] = 0;
                    }
                    else
                    {
                        int obsLeft = isNeighbourObstacle(cell, LEFT);
                        int obsRight = isNeighbourObstacle(cell, RIGHT);
                        V[i][j] = -V[i + obsLeft - obsRight][j];
                    }
                }
                // Compute u
                if (!skipU(cell))
                {
                    if (isNeighbourFluid(cell, RIGHT))
                    {
                        U[i][j] = 0;
                    }
                    else
                    {
                        int obsBottom = isNeighbourObstacle(cell, BOT);
                        int obsTop = isNeighbourObstacle(cell, TOP);
                        U[i][j] = -U[i][j + obsBottom - obsTop];
                    }
                }
            }
            else // if (isFluid(cell))
            {
                //compute V
                if (isNeighbourObstacle(cell, TOP) && (j != jmax))
                {
                    V[i][j] = 0;
                }
                //compute U
                if (isNeighbourObstacle(cell, RIGHT) && (i != imax))
                {
                    U[i][j] = 0;
                }
            }
        }
    }
}

void leftboundary(int imax, int jmax, double **U, double **V, BoundaryInfo bI[4])
{
    for (int j = 1; j <= jmax; j++)
    {
        
        //bI[2] == LEFT
        if (bI[2].typeU == DIRICHLET)
        {
            if (bI[2].constU)
            {
                U[0][j] = *(bI[2].valuesU);
            }
            else
            {
                U[0][j] = (bI[2].valuesU)[j - 1];
            }
        }
        else
        {
            U[0][j] = U[1][j];
        }
        
        if (bI[2].typeV == DIRICHLET)
        {
            if (bI[2].constV)
            {
                V[0][j] = 2 * (bI[2].valuesV)[0] - V[1][j];
            }
            else
            {
                V[0][j] = 2 * (bI[2].valuesV)[j - 1] - V[1][j];
            }
        }
        else
        {
            V[0][j] = V[1][j];
        }
    }
}

void rightboundary(int imax, int jmax, double **U, double **V, BoundaryInfo bI[4])
{
    for (int j = 1; j <= jmax; j++)
    {
        
        //bI[3] == RIGHT
        if (bI[3].typeU == DIRICHLET)
        {
            if (bI[3].constU)
            {
                U[imax][j] = *(bI[3].valuesU);
            }
            else
            {
                U[imax][j] = (bI[3].valuesU)[j - 1];
            }
        }
        else
        {
            U[imax][j] = U[imax - 1][j];
        }
        
        if (bI[3].typeV == DIRICHLET)
        {
            if (bI[3].constV)
            {
                V[imax + 1][j] = 2 * (bI[3].valuesV)[0] - V[imax][j];
            }
            else
            {
                V[imax + 1][j] = 2 * (bI[3].valuesV)[j - 1] - V[imax][j];
            }
        }
        else
        {
            V[imax + 1][j] = V[imax][j];
        }
    }
}

void topboundary(int imax, int jmax, double **U, double **V, BoundaryInfo bI[4])
{
    for (int i = 1; i <= imax; i++)
    {
        
        //bI[0] == TOP
        if (bI[0].typeV == DIRICHLET)
        {
            if (bI[0].constV)
            {
                V[i][jmax] = *(bI[0].valuesV);
            }
            else
            {
                V[i][jmax] = (bI[0].valuesV)[i - 1];
            }
        }
        else
        {
            V[i][jmax] = V[i][jmax - 1];
        }
        
        if (bI[0].typeU == DIRICHLET)
        {
            if (bI[0].constU)
            {
                U[i][jmax + 1] = 2 * (bI[0].valuesU)[0] - U[i][jmax];
            }
            else
            {
                U[i][jmax + 1] = 2 * (bI[0].valuesU)[i - 1] - U[i][jmax];
            }
        }
        else
        {
            U[i][jmax + 1] = U[i][jmax];
        }
    }
}

void bottomboundary(int imax, int jmax, double **U, double **V, BoundaryInfo bI[4])
{
    for (int i = 1; i <= imax; i++)
    {
        
        //bI[1] == BOTTOM
        if (bI[1].typeV == DIRICHLET)
        {
            if (bI[1].constV)
            {
                V[i][0] = *(bI[1].valuesV);
            }
            else
            {
                V[i][0] = (bI[1].valuesV)[i - 1];
            }
        }
        else
        {
            V[i][0] = V[i][1];
        }
        
        if (bI[0].typeU == DIRICHLET)
        {
            if (bI[1].constU)
            {
                U[i][0] = 2 * (bI[1].valuesU)[0] - U[i][1];
            }
            else
            {
                U[i][0] = 2 * (bI[1].valuesU)[i - 1] - U[i][1];
            }
        }
        else
        {
            U[i][0] = U[i][1];
        }
    }
}

void
initBoundaryInfo(BoundaryInfo *boundaryInfo, BoundaryType typeU, BoundaryType typeV, int numValuesU, int numValuesV)
{
    boundaryInfo->typeU = typeU;
    boundaryInfo->typeV = typeV;
    if (typeU == NEUMANN)
    {
        boundaryInfo->constU = 1;
        boundaryInfo->valuesU = NULL;
    }
    else
    {
        boundaryInfo->valuesU = calloc(numValuesU, sizeof(double));
        boundaryInfo->constU = (numValuesU == 1);
    }
    if (typeV == NEUMANN)
    {
        boundaryInfo->constV = 1;
        boundaryInfo->valuesV = NULL;
    }
    else
    {
        boundaryInfo->valuesV = calloc(numValuesV, sizeof(double));
        boundaryInfo->constV = (numValuesV == 1);
    }
}

//eof
