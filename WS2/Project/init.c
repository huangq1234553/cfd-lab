#include "helper.h"
#include "init.h"
#include "boundary_val.h"
#include "logger.h"

int read_parameters(const char *szFileName, double *Re, double *UI, double *VI, double *PI, double *GX, double *GY,
                    double *t_end, double *xlength, double *ylength, double *dt, double *dx, double *dy, int *imax,
                    int *jmax, double *alpha, double *omg, double *tau, int *itermax, double *eps, double *dt_value,
                    char *problem, char *geometry, BoundaryInfo boundaryInfo[4],
                    double *beta, double *TI, double *T_h, double *T_c, double* Pr)    /* path/filename to geometry file */
{
    READ_DOUBLE(szFileName, *xlength, REQUIRED);
    READ_DOUBLE(szFileName, *ylength, REQUIRED);
    
    READ_DOUBLE(szFileName, *Re, REQUIRED);
    READ_DOUBLE(szFileName, *t_end, REQUIRED);
    READ_DOUBLE(szFileName, *dt, REQUIRED);
    
    READ_INT   (szFileName, *imax, REQUIRED);
    READ_INT   (szFileName, *jmax, REQUIRED);
    
    READ_DOUBLE(szFileName, *omg, REQUIRED);
    READ_DOUBLE(szFileName, *eps, REQUIRED);
    READ_DOUBLE(szFileName, *tau, REQUIRED);
    READ_DOUBLE(szFileName, *alpha, REQUIRED);
    
    READ_INT   (szFileName, *itermax, REQUIRED);
    READ_DOUBLE(szFileName, *dt_value, REQUIRED);
    
    READ_DOUBLE(szFileName, *UI, REQUIRED);
    READ_DOUBLE(szFileName, *VI, REQUIRED);
    READ_DOUBLE(szFileName, *GX, REQUIRED);
    READ_DOUBLE(szFileName, *GY, REQUIRED);
    READ_DOUBLE(szFileName, *PI, REQUIRED);
    
    READ_DOUBLE(szFileName, *beta, OPTIONAL);
    READ_DOUBLE(szFileName, *TI, OPTIONAL);
    READ_DOUBLE(szFileName, *T_h, OPTIONAL);
    READ_DOUBLE(szFileName, *T_c, OPTIONAL);
    READ_DOUBLE(szFileName, *Pr, OPTIONAL);

    READ_STRING(szFileName, problem, REQUIRED);
    READ_STRING(szFileName, geometry, REQUIRED);
    
    *dx = *xlength / (double) (*imax);
    *dy = *ylength / (double) (*jmax);
    
    // Now read boundary-related variables
//    initBoundaryInfo(boundaryInfo+LEFTBOUNDARY,NEUMANN,NEUMANN,1,1);
    initBoundaryInfo(boundaryInfo+LEFTBOUNDARY,DIRICHLET,DIRICHLET,1,1);
    *(boundaryInfo[LEFTBOUNDARY].valuesU) = 1;
    *(boundaryInfo[LEFTBOUNDARY].valuesV) = 0;
    //
    initBoundaryInfo(boundaryInfo+RIGHTBOUNDARY,NEUMANN,NEUMANN,1,1); // outflow
//    initBoundaryInfo(boundaryInfo+RIGHTBOUNDARY,DIRICHLET,DIRICHLET,1,1); // no-slip
//    *(boundaryInfo[RIGHTBOUNDARY].valuesU) = 0;
//    *(boundaryInfo[RIGHTBOUNDARY].valuesV) = 0;
    //
    
    initBoundaryInfo(boundaryInfo+TOPBOUNDARY,NEUMANN,NEUMANN,1,1); // outflow
//    initBoundaryInfo(boundaryInfo+TOPBOUNDARY,DIRICHLET,DIRICHLET,1,1);
//    *(boundaryInfo[TOPBOUNDARY].valuesU) = 0;
//    *(boundaryInfo[TOPBOUNDARY].valuesV) = 0;
    initBoundaryInfo(boundaryInfo+BOTTOMBOUNDARY,DIRICHLET,DIRICHLET,1,1);
    *(boundaryInfo[BOTTOMBOUNDARY].valuesU) = 0;
    *(boundaryInfo[BOTTOMBOUNDARY].valuesV) = 0;
    //
    return 1;
}

void init_uvpt(double UI, double VI, double PI, double TI, int imax, int jmax, double **U, double **V, double **P,
               double **T, int **Flags)
{
    init_matrix(U, 0, imax + 1, 0, jmax + 1, UI);
    init_matrix(V, 0, imax + 1, 0, jmax + 1, VI);
    init_matrix(P, 0, imax + 1, 0, jmax + 1, PI);
    init_matrix(T, 0, imax + 1, 0, jmax + 1, TI);
    for (int i = 0; i <= imax+1; ++i)
    {
        for (int j = 0; j <= jmax+1; ++j)
        {
            if (isObstacle(Flags[i][j]))
            {
                U[i][j] = 0;
                V[i][j] = 0;
                P[i][j] = 0;
                T[i][j] = 0;
            }
        }
    }
}

void init_flag(
        char *problem,
        char *geometry,
        int imax,
        int jmax,
        int **Flag,
        int *counter
)
{
    int **pic = NULL;
    
    pic = read_pgm(geometry); // NOTE: this is covering just the inner part of the image, so it is imax*jmax

    // Set the outer boundary + the first inner layers
    for (int i = 0; i < imax + 1; ++i)
    {
        // Outer boundary
        Flag[i][0] = 1;
        Flag[i][jmax + 1] = 1;
    }
    for (int j = 0; j < jmax + 1; ++j)
    {
        // Outer boundary
        Flag[0][j] = 1;
        Flag[imax + 1][j] = 1;
    }
    // Set the inner domain geometry
    for (int i = 1; i < imax+1; i++)
    {
        for (int j = 1; j < jmax+1; j++)
        {
            Flag[i][j] = pic[i - 1][j - 1];
        }
    }
    *counter = 0;
    // Set the boundary domain flags
    for (int j=1; j<jmax+1; ++j)
    {
        Flag[0][j] += (1<<TOP) * 1
                      + (1<<BOT) * 1
                      + (1<<LEFT) * 1
                      + (1<<RIGHT) * isObstacle(Flag[1][j]);
        Flag[imax+1][j] += (1<<TOP) * 1
                      + (1<<BOT) * 1
                      + (1<<LEFT) * isObstacle(Flag[imax][j])
                      + (1<<RIGHT) * 1;
    }
    for (int i=1; i<imax+1; ++i)
    {
        Flag[i][0] += (1<<TOP) * isObstacle(Flag[i][1])
                      + (1<<BOT) * 1
                      + (1<<LEFT) * 1
                      + (1<<RIGHT) * 1;
        Flag[i][jmax+1] += (1<<TOP) * 1
                           + (1<<BOT) * isObstacle(Flag[i][jmax])
                           + (1<<LEFT) * 1
                           + (1<<RIGHT) * 1;
    }
    // Set the inner domain flags
    for (int j = jmax; j > 0; j--)
    {
        for (int i = 1; i < imax + 1; i++)
        {
            Flag[i][j] += (1<<TOP) * isObstacle(Flag[i][j + 1])
                          + (1<<BOT) * isObstacle(Flag[i][j - 1])
                          + (1<<LEFT) * isObstacle(Flag[i - 1][j])
                          + (1<<RIGHT) * isObstacle(Flag[i + 1][j]);
            logRawString("%d ", isObstacle(Flag[i][j]));
            (*counter) += isFluid(Flag[i][j]);
        }
        logRawString("\n");
    }
    logMsg("Total fluid cells in domain: %d", (*counter));
    geometryCheck(Flag, imax, jmax);
    free_imatrix(pic, 0, imax + 1, 0, jmax + 1);
}
