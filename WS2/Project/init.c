#include "helper.h"
#include "init.h"
#include "logger.h"
#include "boundary_configurator.h"

void setDefaultStringIfRequired(char *variable, const char *defaultValue)
{
    if (strcmp(variable, "NULLSTRING") == 0)
    {
        strcpy(variable, defaultValue);
    }
}

int read_parameters(const char *szFileName, double *Re, double *UI, double *VI, double *PI, double *GX, double *GY,
                    double *t_end, double *xlength, double *ylength, double *dt, double *dx, double *dy, int *imax,
                    int *jmax, double *alpha, double *omg, double *tau, int *itermax, double *eps, double *dt_value,
                    char *problem, char *geometry, BoundaryInfo boundaryInfo[4],
                    double *beta, double *TI, double *T_h, double *T_c,
                    double *Pr)    /* path/filename to geometry file */
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
//    READ_DOUBLE(szFileName, *T_h, OPTIONAL);  // NOT REQUIRED
//    READ_DOUBLE(szFileName, *T_c, OPTIONAL);  // NOT REQUIRED
    READ_DOUBLE(szFileName, *Pr, OPTIONAL);
    
    READ_STRING(szFileName, problem, REQUIRED);
    READ_STRING(szFileName, geometry, REQUIRED);
    
    *dx = *xlength / (double) (*imax);
    *dy = *ylength / (double) (*jmax);
    
    // Now read boundary-related variables
    char left_boundary_type[16];
    char left_boundary_temp_type[16];
    char right_boundary_type[16];
    char right_boundary_temp_type[16];
    char top_boundary_type[16];
    char top_boundary_temp_type[16];
    char bottom_boundary_type[16];
    char bottom_boundary_temp_type[16];
    double left_boundary_U;
    double left_boundary_V;
    double left_boundary_T;
    double left_boundary_qN;
    double left_boundary_k;
    double right_boundary_U;
    double right_boundary_V;
    double right_boundary_T;
    double right_boundary_qN;
    double right_boundary_k;
    double top_boundary_U;
    double top_boundary_V;
    double top_boundary_T;
    double top_boundary_qN;
    double top_boundary_k;
    double bottom_boundary_U;
    double bottom_boundary_V;
    double bottom_boundary_T;
    double bottom_boundary_qN;
    double bottom_boundary_k;
    
    char *boundaryTypeDefault = "NOSLIP";
    char *boundaryTempTypeDefault = "NEUMANN";
    
    READ_STRING(szFileName, left_boundary_type, OPTIONAL);
    setDefaultStringIfRequired(left_boundary_type, boundaryTypeDefault);
    READ_STRING(szFileName, left_boundary_temp_type, OPTIONAL);
    setDefaultStringIfRequired(left_boundary_temp_type, boundaryTempTypeDefault);
    
    READ_STRING(szFileName, right_boundary_type, OPTIONAL);
    setDefaultStringIfRequired(right_boundary_type, boundaryTypeDefault);
    READ_STRING(szFileName, right_boundary_temp_type, OPTIONAL);
    setDefaultStringIfRequired(right_boundary_temp_type, boundaryTempTypeDefault);
    
    READ_STRING(szFileName, top_boundary_type, OPTIONAL);
    setDefaultStringIfRequired(top_boundary_type, boundaryTypeDefault);
    READ_STRING(szFileName, top_boundary_temp_type, OPTIONAL);
    setDefaultStringIfRequired(top_boundary_temp_type, boundaryTempTypeDefault);
    
    READ_STRING(szFileName, bottom_boundary_type, OPTIONAL);
    setDefaultStringIfRequired(bottom_boundary_type, boundaryTypeDefault);
    READ_STRING(szFileName, bottom_boundary_temp_type, OPTIONAL);
    setDefaultStringIfRequired(bottom_boundary_temp_type, boundaryTempTypeDefault);
        
    READ_DOUBLE(szFileName, left_boundary_U, OPTIONAL);
    READ_DOUBLE(szFileName, left_boundary_V, OPTIONAL);
    READ_DOUBLE(szFileName, left_boundary_T, OPTIONAL);
    READ_DOUBLE(szFileName, left_boundary_qN, OPTIONAL);
    READ_DOUBLE(szFileName, left_boundary_k, OPTIONAL);
    if (left_boundary_k == 0.0)
        left_boundary_k = 1;
    
    READ_DOUBLE(szFileName, right_boundary_U, OPTIONAL);
    READ_DOUBLE(szFileName, right_boundary_V, OPTIONAL);
    READ_DOUBLE(szFileName, right_boundary_T, OPTIONAL);
    READ_DOUBLE(szFileName, right_boundary_qN, OPTIONAL);
    READ_DOUBLE(szFileName, right_boundary_k, OPTIONAL);
    if (right_boundary_k == 0.0)
        right_boundary_k = 1;
    
    READ_DOUBLE(szFileName, top_boundary_U, OPTIONAL);
    READ_DOUBLE(szFileName, top_boundary_V, OPTIONAL);
    READ_DOUBLE(szFileName, top_boundary_T, OPTIONAL);
    READ_DOUBLE(szFileName, top_boundary_qN, OPTIONAL);
    READ_DOUBLE(szFileName, top_boundary_k, OPTIONAL);
    if (top_boundary_k == 0.0)
        top_boundary_k = 1;
    
    READ_DOUBLE(szFileName, bottom_boundary_U, OPTIONAL);
    READ_DOUBLE(szFileName, bottom_boundary_V, OPTIONAL);
    READ_DOUBLE(szFileName, bottom_boundary_T, OPTIONAL);
    READ_DOUBLE(szFileName, bottom_boundary_qN, OPTIONAL);
    READ_DOUBLE(szFileName, bottom_boundary_k, OPTIONAL);
    if (bottom_boundary_k == 0.0)
        bottom_boundary_k = 1;
    
    configureBoundary(boundaryInfo, LEFTBOUNDARY, left_boundary_type, left_boundary_temp_type, left_boundary_U, left_boundary_V,
                      left_boundary_T, left_boundary_qN, left_boundary_k, *dx);
    configureBoundary(boundaryInfo, RIGHTBOUNDARY, right_boundary_type, right_boundary_temp_type, right_boundary_U, right_boundary_V,
                      right_boundary_T, right_boundary_qN, right_boundary_k, *dy);
    configureBoundary(boundaryInfo, TOPBOUNDARY, top_boundary_type, top_boundary_temp_type, top_boundary_U, top_boundary_V,
                      top_boundary_T, top_boundary_qN, top_boundary_k, *dx);
    configureBoundary(boundaryInfo, BOTTOMBOUNDARY, bottom_boundary_type, bottom_boundary_temp_type, bottom_boundary_U, bottom_boundary_V,
                      bottom_boundary_T, bottom_boundary_qN, bottom_boundary_k, *dy);
   
    // TODO: add support for more complex profiles and/or autogeneration of parabolic one. Do this into the new boundary_configurator.c file
    
    return 1;
}

void init_uvpt(double UI, double VI, double PI, double TI, int imax, int jmax, double **U, double **V, double **P,
               double **T, int **Flags)
{
    init_matrix(U, 0, imax + 1, 0, jmax + 1, UI);
    init_matrix(V, 0, imax + 1, 0, jmax + 1, VI);
    init_matrix(P, 0, imax + 1, 0, jmax + 1, PI);
    init_matrix(T, 0, imax + 1, 0, jmax + 1, TI);
    for (int i = 1; i <= imax; ++i)
    {
        for (int j = 1; j <= jmax; ++j)
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
    for (int i = 1; i < imax + 1; i++)
    {
        for (int j = 1; j < jmax + 1; j++)
        {
            Flag[i][j] = pic[i - 1][j - 1];
        }
    }
    *counter = 0;
    // Set the boundary domain flags
    for (int j = 1; j < jmax + 1; ++j)
    {
        Flag[0][j] += (1 << TOP) * 1
                      + (1 << BOT) * 1
                      + (1 << LEFT) * 1
                      + (1 << RIGHT) * isObstacle(Flag[1][j]);
        Flag[imax + 1][j] += (1 << TOP) * 1
                             + (1 << BOT) * 1
                             + (1 << LEFT) * isObstacle(Flag[imax][j])
                             + (1 << RIGHT) * 1;
    }
    for (int i = 1; i < imax + 1; ++i)
    {
        Flag[i][0] += (1 << TOP) * isObstacle(Flag[i][1])
                      + (1 << BOT) * 1
                      + (1 << LEFT) * 1
                      + (1 << RIGHT) * 1;
        Flag[i][jmax + 1] += (1 << TOP) * 1
                             + (1 << BOT) * isObstacle(Flag[i][jmax])
                             + (1 << LEFT) * 1
                             + (1 << RIGHT) * 1;
    }
    // Set the inner domain flags
    for (int j = jmax; j > 0; j--)
    {
        for (int i = 1; i < imax + 1; i++)
        {
            Flag[i][j] += (1 << TOP) * isObstacle(Flag[i][j + 1])
                          + (1 << BOT) * isObstacle(Flag[i][j - 1])
                          + (1 << LEFT) * isObstacle(Flag[i - 1][j])
                          + (1 << RIGHT) * isObstacle(Flag[i + 1][j]);
            (*counter) += isFluid(Flag[i][j]);
        }
    }
    logMsg("Total fluid cells in domain: %d", (*counter));
    geometryCheck(Flag, imax, jmax);
    free_imatrix(pic, 0, imax + 1, 0, jmax + 1);
}
