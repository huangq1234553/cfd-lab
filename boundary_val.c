#include "boundary_val.h"
#include "helper.h"
#include "logger.h"

void boundaryval(int imax, int jmax, double **U, double **V, double **T, int **Flags,
                 BoundaryInfo *boundaryInfo)
{
    // Setting boundary conditions on the outer boundary
    setLeftBoundaryValues(imax, jmax, U, V, T, Flags, boundaryInfo);
    setRightBoundaryValues(imax, jmax, U, V, T, Flags, boundaryInfo);
    setTopBoundaryValues(imax, jmax, U, V, T, Flags, boundaryInfo);
    setBottomBoundaryValues(imax, jmax, U, V, T, Flags, boundaryInfo);
    
    // Boundary values at geometries in the internal part of the domain
    for (int i = 1; i <= imax; ++i)
    {
        for (int j = 1; j <= jmax; ++j)
        {
            setEdgeBoundaryValues(imax, jmax, U, V, Flags, i, j);
            setCenterBoundaryValues(imax, jmax, T, Flags, i, j, true);
        }
    }
}

void setLeftBoundaryValues(int imax, int jmax, double **U, double **V, double **T, int **Flags,
                           BoundaryInfo *boundaryInfo)
{
    int Flag;
    for (int j = 1; j <= jmax; j++)
    {
        // Set the velocity boundary values
        Flag = Flags[0][j];
        int rightNeighbourIsObstacle = isNeighbourObstacle(Flags[0][j],RIGHT);
        if(rightNeighbourIsObstacle){
            ;
        }
        else {
            if (Flag >> NSBIT & 1 || Flag >> CBIT & 1) {
                U[0][j] = 0;
                V[0][j] = -V[1][j];
            } else if (Flag >> FSBIT & 1) {
                U[0][j] = 0;
                V[0][j] = V[1][j];
            } else if (Flag >> IFBIT & 1) {
                U[0][j] = (boundaryInfo[LEFTBOUNDARY].valuesDirichletU)[0];
                V[0][j] = 2 * (boundaryInfo[LEFTBOUNDARY].valuesDirichletV)[0] - V[1][j];
            } else {
                U[0][j] = U[1][j];
                V[0][j] = V[1][j];
            }
            // Set temperature boundary values
            if (Flag >> TBIT & 1) {
                T[0][j] = 2 * (boundaryInfo[LEFTBOUNDARY].valuesDirichletT)[0] - T[1][j];
            } else {
                T[0][j] = T[1][j] + boundaryInfo[LEFTBOUNDARY].coeff;
            }
        }
    }
}

void setRightBoundaryValues(int imax, int jmax, double **U, double **V, double **T, int **Flags,
                            BoundaryInfo *boundaryInfo) {
    int Flag;
    for (int j = 1; j <= jmax; j++) {
        Flag = Flags[imax + 1][j];
        int leftNeighbourIsObstacle = isNeighbourObstacle(Flags[imax + 1][j], LEFT);
        if (leftNeighbourIsObstacle) { ;
        }
        else {
            if (Flag >> NSBIT & 1 || Flag >> CBIT & 1) {
                U[imax][j] = 0;
                V[imax + 1][j] = -V[imax][j];
            } else if (Flag >> FSBIT & 1) {
                U[imax][j] = 0;
                V[imax + 1][j] = V[imax][j];
            } else if (Flag >> IFBIT & 1) {
                U[imax][j] = (boundaryInfo[RIGHTBOUNDARY].valuesDirichletU)[0];
                V[imax][j] = 2 * (boundaryInfo[RIGHTBOUNDARY].valuesDirichletV)[0] - V[imax][j];
            } else {
                U[imax][j] = U[imax - 1][j];
                V[imax + 1][j] = V[imax][j];
            }

            if (Flag >> TBIT & 1) {
                T[imax + 1][j] = 2 * (boundaryInfo[RIGHTBOUNDARY].valuesDirichletT)[0] - T[imax][j];
            } else {
                T[imax + 1][j] = T[imax][j] + boundaryInfo[RIGHTBOUNDARY].coeff;
            }
        }
    }
}

void setTopBoundaryValues(int imax, int jmax, double **U, double **V, double **T, int **Flags,
                          BoundaryInfo *boundaryInfo)
{
    int Flag;
    for (int i = 1; i <= imax; i++)
    {
        Flag = Flags[i][jmax+1];
        int bottomNeighbourIsObstacle = isNeighbourObstacle(Flags[i][jmax+1],BOT);
        if (bottomNeighbourIsObstacle) { ;
        }
        else {
            if (Flag >> NSBIT & 1 || Flag >> CBIT & 1) {
                U[i][jmax+1] = -U[i][jmax];
                V[i][jmax] = 0;
            } else if (Flag >> FSBIT & 1) {
                U[i][jmax+1] = U[i][jmax];
                V[i][jmax] = 0;
            } else if (Flag >> IFBIT & 1) {
                U[i][jmax+1] = 2 * (boundaryInfo[TOPBOUNDARY].valuesDirichletU)[0] - U[i][jmax];
                V[i][jmax] = (boundaryInfo[TOPBOUNDARY].valuesDirichletV)[0] ;
            } else {
                U[i][jmax+1] = U[i][jmax];
                V[i][jmax] = V[imax][jmax-1];
            }

            if (Flag >> TBIT & 1) {
                T[i][jmax+1] = 2 * (boundaryInfo[TOPBOUNDARY].valuesDirichletT)[0] - T[i][jmax];
            } else {
                T[i][jmax+1] = T[i][jmax] + boundaryInfo[TOPBOUNDARY].coeff;
            }
        }
    }
}

void setBottomBoundaryValues(int imax, int jmax, double **U, double **V, double **T, int **Flags,
                             BoundaryInfo *boundaryInfo)
{
    int Flag;
    for (int i = 1; i <= imax; i++)
    {
        Flag = Flags[i][0];
        int topNeighbourIsObstacle = isNeighbourObstacle(Flags[i][0],TOP);
        if (topNeighbourIsObstacle) { ;
        }
        else {
            if (Flag >> NSBIT & 1 || Flag >> CBIT & 1) {
                U[i][0] = -U[i][1];
                V[i][0] = 0;
            } else if (Flag >> FSBIT & 1) {
                U[i][0] = U[i][1];
                V[i][0] = 0;
            } else if (Flag >> IFBIT & 1) {
                U[i][0] = 2 * (boundaryInfo[TOPBOUNDARY].valuesDirichletU)[0] - U[i][1];
                V[i][0] = (boundaryInfo[TOPBOUNDARY].valuesDirichletV)[0] ;
            } else {
                U[i][0] = U[i][1];
                V[i][0] = V[i][1];
            }

            if (Flag >> TBIT & 1) {
                T[i][0] = 2 * (boundaryInfo[TOPBOUNDARY].valuesDirichletT)[0] - T[i][1];
            } else {
                T[i][0] = T[i][0] + boundaryInfo[BOTTOMBOUNDARY].coeff;
            }
        }
    }
}

/*void initBoundaryInfo(BoundaryInfo *boundaryInfo, HLBoundaryType hlBoundaryType, BoundaryType typeU, BoundaryType typeV,
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
}*/

void setEdgeBoundaryValues(int imax, int jmax, double **U, double **V, int  **Flags, int i, int j)
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

void setCenterBoundaryValues(int imax, int jmax, double **Q, int **Flags, int i, int j, bool isCoupled)
{
    int C = Flags[i][j];
    // proceed if obstacle
    if (isObstacle(C) && !(isCoupling(C) && isCoupled))
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

void freeAllBoundaryInfo(BoundaryInfo boundaryInfo[4])
{
    for (int i=0; i<4; ++i)
    {
        free(boundaryInfo[i].valuesDirichletU);
        free(boundaryInfo[i].valuesDirichletV);
        free(boundaryInfo[i].valuesDirichletT);
    }
}

//eof
