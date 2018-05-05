#include "helper.h"
#include "init.h"
#include "boundary_val.h"

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
    initBoundaryInfo(boundaryInfo+L,DIRICHLET,DIRICHLET,1,1);
    *(boundaryInfo[L].valuesU) = 0;
    *(boundaryInfo[L].valuesV) = 0;
    initBoundaryInfo(boundaryInfo+R,DIRICHLET,DIRICHLET,1,1);
    *(boundaryInfo[R].valuesU) = 0;
    *(boundaryInfo[R].valuesV) = 0;
    initBoundaryInfo(boundaryInfo+T,DIRICHLET,DIRICHLET,1,1);
    *(boundaryInfo[T].valuesU) = 1;
    *(boundaryInfo[T].valuesV) = 0;
    initBoundaryInfo(boundaryInfo+B,DIRICHLET,DIRICHLET,1,1);
    *(boundaryInfo[B].valuesU) = 0;
    *(boundaryInfo[B].valuesV) = 0;
    //
    return 1;
}

void init_uvpt(
        double UI,
        double VI,
        double PI,
        double TI,
        int imax,
        int jmax,
        double **U,
        double **V,
        double **P,
        double **T
)
{
    init_matrix(U, 0, imax + 1, 0, jmax + 1, UI);
    init_matrix(V, 0, imax + 1, 0, jmax + 1, VI);
    init_matrix(P, 0, imax + 1, 0, jmax + 1, PI);
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
    for (int j = 0; j < imax + 1; ++j)
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
    // Set the inner domain flags
    for (int j = jmax; j > 0; j--)
    {
        for (int i = 1; i < imax + 1; i++)
        {
            Flag[i][j] += (1<<TOP) * isObstacle(Flag[i][j + 1])
                          + (1<<BOT) * isObstacle(Flag[i][j - 1])
                          + (1<<LEFT) * isObstacle(Flag[i - 1][j])
                          + (1<<RIGHT) * isObstacle(Flag[i + 1][j]);
            printf("%d ", isCorner(Flag[i][j]));
            (*counter) += isFluid(Flag[i][j]);
        }
        printf("\n");
    }
    printf("%d\n", (*counter));
    geometryCheck(Flag, imax, jmax);
    free_imatrix(pic, 0, imax + 1, 0, jmax + 1);
}
