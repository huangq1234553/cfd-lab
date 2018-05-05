#include "boundary_val.h"
#include "helper.h"

void boundaryvalues(int imax, int jmax, double **U, double **V, int **Flags, BoundaryInfo boundaryInfo[4])
{
    // Setting boundary conditions on the outer boundary
    setLeftBoundaryVelocities(imax, jmax, U, V, boundaryInfo);
    setRightBoundaryVelocities(imax, jmax, U, V, boundaryInfo);
    setTopBoundaryVelocities(imax, jmax, U, V, boundaryInfo);
    setBottomBoundaryVelocities(imax, jmax, U, V, boundaryInfo);
    
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

void setLeftBoundaryVelocities(int imax, int jmax, double **U, double **V, BoundaryInfo *bI)
{
    for (int j = 1; j <= jmax; j++)
    {
        //bI[2] == LEFT
        if (bI[L].typeU == DIRICHLET)
        {
            if (bI[L].constU)
            {
                U[0][j] = *(bI[L].valuesU);
            }
            else
            {
                U[0][j] = (bI[L].valuesU)[j - 1];
            }
        }
        else
        {
            U[0][j] = U[1][j];
        }
        
        if (bI[L].typeV == DIRICHLET)
        {
            if (bI[L].constV)
            {
                V[0][j] = 2 * (bI[L].valuesV)[0] - V[1][j];
            }
            else
            {
                V[0][j] = 2 * (bI[L].valuesV)[j - 1] - V[1][j];
            }
        }
        else
        {
            V[0][j] = V[1][j];
        }
    }
}

void setRightBoundaryVelocities(int imax, int jmax, double **U, double **V, BoundaryInfo *bI)
{
    for (int j = 1; j <= jmax; j++)
    {
        //bI[3] == RIGHT
        if (bI[R].typeU == DIRICHLET)
        {
            if (bI[R].constU)
            {
                U[imax][j] = *(bI[R].valuesU);
            }
            else
            {
                U[imax][j] = (bI[R].valuesU)[j - 1];
            }
        }
        else
        {
            U[imax][j] = U[imax - 1][j];
        }
        
        if (bI[R].typeV == DIRICHLET)
        {
            if (bI[R].constV)
            {
                V[imax + 1][j] = 2 * (bI[R].valuesV)[0] - V[imax][j];
            }
            else
            {
                V[imax + 1][j] = 2 * (bI[R].valuesV)[j - 1] - V[imax][j];
            }
        }
        else
        {
            V[imax + 1][j] = V[imax][j];
        }
    }
}

void setTopBoundaryVelocities(int imax, int jmax, double **U, double **V, BoundaryInfo *bI)
{
    for (int i = 1; i <= imax; i++)
    {
        
        //bI[0] == TOP
        if (bI[T].typeV == DIRICHLET)
        {
            if (bI[T].constV)
            {
                V[i][jmax] = *(bI[T].valuesV);
            }
            else
            {
                V[i][jmax] = (bI[T].valuesV)[i - 1];
            }
        }
        else
        {
            V[i][jmax] = V[i][jmax - 1];
        }
        
        if (bI[T].typeU == DIRICHLET)
        {
            if (bI[T].constU)
            {
                U[i][jmax + 1] = 2 * (bI[T].valuesU)[0] - U[i][jmax];
            }
            else
            {
                U[i][jmax + 1] = 2 * (bI[T].valuesU)[i - 1] - U[i][jmax];
            }
        }
        else
        {
            U[i][jmax + 1] = U[i][jmax];
        }
    }
}

void setBottomBoundaryVelocities(int imax, int jmax, double **U, double **V, BoundaryInfo *bI)
{
    for (int i = 1; i <= imax; i++)
    {
        
        //bI[1] == BOTTOM
        if (bI[B].typeV == DIRICHLET)
        {
            if (bI[B].constV)
            {
                V[i][0] = *(bI[B].valuesV);
            }
            else
            {
                V[i][0] = (bI[B].valuesV)[i - 1];
            }
        }
        else
        {
            V[i][0] = V[i][1];
        }
        
        if (bI[B].typeU == DIRICHLET)
        {
            if (bI[B].constU)
            {
                U[i][0] = 2 * (bI[B].valuesU)[0] - U[i][1];
            }
            else
            {
                U[i][0] = 2 * (bI[B].valuesU)[i - 1] - U[i][1];
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
