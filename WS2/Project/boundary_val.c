#include "boundary_val.h"
#include "helper.h"
#include "logger.h"

void boundaryval(int imax, int jmax, double **U, double **V, double **T, int **Flags,
                 BoundaryInfo *boundaryInfo)
{
    // Setting boundary conditions on the outer boundary
    setLeftBoundaryVelocities(imax, jmax, U, V, T, Flags, boundaryInfo);
    setRightBoundaryVelocities(imax, jmax, U, V, T, Flags, boundaryInfo);
    setTopBoundaryVelocities(imax, jmax, U, V, T, Flags, boundaryInfo);
    setBottomBoundaryVelocities(imax, jmax, U, V, T, Flags, boundaryInfo);
    
    // Boundary values at geometries in the internal part of the domain
    for (int i = 1; i <= imax; ++i)
    {
        for (int j = 1; j <= jmax; ++j)
        {
            setEdgeBoundaryValues(imax, jmax, U, V, Flags, i, j);
            setCenterBoundaryValues(imax, jmax, T, Flags, i, j);
        }
    }
}

void setLeftBoundaryVelocities(int imax, int jmax, double **U, double **V, double **T, int **Flags,
                               BoundaryInfo *boundaryInfo)
{
    for (int j = 1; j <= jmax; j++)
    {
        // Set the velocity boundary values
        int rightNeighbourIsFluid = isNeighbourFluid(Flags[0][j],RIGHT);
        // U values
        if (boundaryInfo[LEFTBOUNDARY].typeU == DIRICHLET
            && rightNeighbourIsFluid)
        {
            if (boundaryInfo[LEFTBOUNDARY].constU)
            {
                U[0][j] = *(boundaryInfo[LEFTBOUNDARY].valuesU);
            }
            else
            {
                U[0][j] = (boundaryInfo[LEFTBOUNDARY].valuesU)[j - 1];
            }
        }
        else
        {
            U[0][j] = U[1][j];
        }
        // V values
        if (boundaryInfo[LEFTBOUNDARY].typeV == DIRICHLET
            && rightNeighbourIsFluid)
        {
            if (boundaryInfo[LEFTBOUNDARY].constV)
            {
                V[0][j] = 2 * (boundaryInfo[LEFTBOUNDARY].valuesV)[0] - V[1][j];
            }
            else
            {
                V[0][j] = 2 * (boundaryInfo[LEFTBOUNDARY].valuesV)[j - 1] - V[1][j];
            }
        }
        else
        {
            V[0][j] = V[1][j];
        }
        // Set temperature boundary values
        if (boundaryInfo[LEFTBOUNDARY].typeT == DIRICHLET
                && rightNeighbourIsFluid)
        {
            if (boundaryInfo[LEFTBOUNDARY].constT)
            {
                T[0][j] = 2 * (boundaryInfo[LEFTBOUNDARY].valuesT)[0] - T[1][j];
            }
            else
            {
                T[0][j] = 2 * (boundaryInfo[LEFTBOUNDARY].valuesT)[j - 1] - T[1][j];
            }
        }
        else
        {
            T[0][j] = T[1][j] + boundaryInfo[LEFTBOUNDARY].coeff;
        }
    }
}

void setRightBoundaryVelocities(int imax, int jmax, double **U, double **V, double **T, int **Flags,
                                BoundaryInfo *boundaryInfo)
{
    for (int j = 1; j <= jmax; j++)
    {
        //boundaryInfo[3] == RIGHT
        int leftNeighbourIsFluid = isNeighbourFluid(Flags[imax+1][j],LEFT);
        if (boundaryInfo[RIGHTBOUNDARY].typeU == DIRICHLET
                && leftNeighbourIsFluid)
        {
            if (boundaryInfo[RIGHTBOUNDARY].constU)
            {
                U[imax][j] = *(boundaryInfo[RIGHTBOUNDARY].valuesU);
            }
            else
            {
                U[imax][j] = (boundaryInfo[RIGHTBOUNDARY].valuesU)[j - 1];
            }
        }
        else
        {
            U[imax][j] = U[imax - 1][j];
        }
        
        if (boundaryInfo[RIGHTBOUNDARY].typeV == DIRICHLET
                && leftNeighbourIsFluid)
        {
            if (boundaryInfo[RIGHTBOUNDARY].constV)
            {
                V[imax + 1][j] = 2 * (boundaryInfo[RIGHTBOUNDARY].valuesV)[0] - V[imax][j];
            }
            else
            {
                V[imax + 1][j] = 2 * (boundaryInfo[RIGHTBOUNDARY].valuesV)[j - 1] - V[imax][j];
            }
        }
        else
        {
            V[imax + 1][j] = V[imax][j];
        }
        // Set temperature boundary values
        if (boundaryInfo[RIGHTBOUNDARY].typeT == DIRICHLET
            && leftNeighbourIsFluid)
        {
            if (boundaryInfo[RIGHTBOUNDARY].constT)
            {
                T[imax+1][j] = 2 * (boundaryInfo[RIGHTBOUNDARY].valuesT)[0] - T[imax][j];
            }
            else
            {
                T[imax+1][j] = 2 * (boundaryInfo[RIGHTBOUNDARY].valuesT)[j - 1] - T[imax][j];
            }
        }
        else
        {
            T[imax+1][j] = T[imax][j] - boundaryInfo[RIGHTBOUNDARY].coeff; // Note: this minus is to keep our convention for coeff sign
        }
    }
}

void setTopBoundaryVelocities(int imax, int jmax, double **U, double **V, double **T, int **Flags,
                              BoundaryInfo *boundaryInfo)
{
    for (int i = 1; i <= imax; i++)
    {
        //boundaryInfo[0] == TOP
        int bottomNeighbourIsFluid = isNeighbourFluid(Flags[i][jmax+1],BOT);
        if (boundaryInfo[TOPBOUNDARY].typeV == DIRICHLET
                && bottomNeighbourIsFluid)
        {
            if (boundaryInfo[TOPBOUNDARY].constV)
            {
                V[i][jmax] = *(boundaryInfo[TOPBOUNDARY].valuesV);
            }
            else
            {
                V[i][jmax] = (boundaryInfo[TOPBOUNDARY].valuesV)[i - 1];
            }
        }
        else
        {
            V[i][jmax] = V[i][jmax - 1];
        }
        
        if (boundaryInfo[TOPBOUNDARY].typeU == DIRICHLET
                && bottomNeighbourIsFluid)
        {
            if (boundaryInfo[TOPBOUNDARY].constU)
            {
                U[i][jmax + 1] = 2 * (boundaryInfo[TOPBOUNDARY].valuesU)[0] - U[i][jmax];
            }
            else
            {
                U[i][jmax + 1] = 2 * (boundaryInfo[TOPBOUNDARY].valuesU)[i - 1] - U[i][jmax];
            }
        }
        else
        {
            U[i][jmax + 1] = U[i][jmax];
        }
        // Set temperature boundary values
        if (boundaryInfo[TOPBOUNDARY].typeT == DIRICHLET
            && bottomNeighbourIsFluid)
        {
            if (boundaryInfo[TOPBOUNDARY].constT)
            {
                T[i][jmax+1] = 2 * (boundaryInfo[TOPBOUNDARY].valuesT)[0] - T[i][jmax];
            }
            else
            {
                T[i][jmax+1] = 2 * (boundaryInfo[TOPBOUNDARY].valuesT)[i - 1] - T[i][jmax];
            }
        }
        else
        {
            T[i][jmax+1] = T[i][jmax] - boundaryInfo[TOPBOUNDARY].coeff; // Note: this minus is to keep our convention for coeff sign
        }
    }
}

void setBottomBoundaryVelocities(int imax, int jmax, double **U, double **V, double **T, int **Flags,
                                 BoundaryInfo *boundaryInfo)
{
    for (int i = 1; i <= imax; i++)
    {
        //boundaryInfo[1] == BOTTOM
        int topNeighbourIsFluid = isNeighbourFluid(Flags[i][0],TOP);
        if (boundaryInfo[BOTTOMBOUNDARY].typeV == DIRICHLET
                && topNeighbourIsFluid)
        {
            if (boundaryInfo[BOTTOMBOUNDARY].constV)
            {
                V[i][0] = *(boundaryInfo[BOTTOMBOUNDARY].valuesV);
            }
            else
            {
                V[i][0] = (boundaryInfo[BOTTOMBOUNDARY].valuesV)[i - 1];
            }
        }
        else
        {
            V[i][0] = V[i][1];
        }
        
        if (boundaryInfo[BOTTOMBOUNDARY].typeU == DIRICHLET
                && topNeighbourIsFluid)
        {
            if (boundaryInfo[BOTTOMBOUNDARY].constU)
            {
                U[i][0] = 2 * (boundaryInfo[BOTTOMBOUNDARY].valuesU)[0] - U[i][1];
            }
            else
            {
                U[i][0] = 2 * (boundaryInfo[BOTTOMBOUNDARY].valuesU)[i - 1] - U[i][1];
            }
        }
        else
        {
            U[i][0] = U[i][1];
        }
        // Set temperature boundary values
        if (boundaryInfo[BOTTOMBOUNDARY].typeT == DIRICHLET
            && topNeighbourIsFluid)
        {
            if (boundaryInfo[BOTTOMBOUNDARY].constT)
            {
                T[i][0] = 2 * (boundaryInfo[BOTTOMBOUNDARY].valuesT)[0] - T[i][1];
            }
            else
            {
                T[i][0] = 2 * (boundaryInfo[BOTTOMBOUNDARY].valuesT)[i - 1] - T[i][1];
            }
        }
        else
        {
            T[i][0] = T[i][1] + boundaryInfo[BOTTOMBOUNDARY].coeff;
        }
    }
}

void initBoundaryInfo(BoundaryInfo *boundaryInfo, BoundaryType typeU, BoundaryType typeV,
                 int numValuesU, int numValuesV)
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
        boundaryInfo->valuesU = calloc((size_t) numValuesU, sizeof(double));
        boundaryInfo->constU = (numValuesU == 1);
    }
    if (typeV == NEUMANN)
    {
        boundaryInfo->constV = 1;
        boundaryInfo->valuesV = NULL;
    }
    else
    {
        boundaryInfo->valuesV = calloc((size_t) numValuesV, sizeof(double));
        boundaryInfo->constV = (numValuesV == 1);
    }
    // Set sane defaults for T-related vars
    boundaryInfo->typeT = NEUMANN;
    boundaryInfo->constT = 1;
    boundaryInfo->valuesT = calloc(1, sizeof(double));
    boundaryInfo->coeff = 0;
}

void setEdgeBoundaryValues(int imax, int jmax, double *const *U, double *const *V, int *const *Flags, int i, int j)
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

void setCenterBoundaryValues(int imax, int jmax, double **Q, int **Flags, int i, int j)
{
    int C = Flags[i][j];
    // proceed if obstacle
    if (isObstacle(C))
    {
        if (isCorner(C))
        {
            Q[i][j] = (Q[i + isNeighbourObstacle(C, LEFT) - isNeighbourObstacle(C, RIGHT)][j] +
                       Q[i][j + isNeighbourObstacle(C, BOT) - isNeighbourObstacle(C, TOP)]) / 2;
        }
        else
        {
            Q[i][j] = (!isNeighbourObstacle(C, TOP)) * Q[i][j + 1];
            Q[i][j] += (!isNeighbourObstacle(C, BOT)) * Q[i][j - 1];
            Q[i][j] += (!isNeighbourObstacle(C, RIGHT)) * Q[i + 1][j];
            Q[i][j] += (!isNeighbourObstacle(C, LEFT)) * Q[i - 1][j];
        }
    }
}

//eof
