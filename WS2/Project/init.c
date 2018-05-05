#include "helper.h"
#include "init.h"

int read_parameters(const char *szFileName,    /* name of the file */
                    double *Re,                 /* reynolds number   */
                    double *UI,                 /* velocity x-direction */
                    double *VI,                 /* velocity y-direction */
                    double *PI,                 /* pressure */
                    double *GX,                 /* gravitation x-direction */
                    double *GY,                 /* gravitation y-direction */
                    double *t_end,              /* end time */
                    double *xlength,            /* length of the domain x-dir.*/
                    double *ylength,            /* length of the domain y-dir.*/
                    double *dt,                 /* time step */
                    double *dx,                 /* length of a cell x-dir. */
                    double *dy,                 /* length of a cell y-dir. */
                    int *imax,                 /* number of cells x-direction*/
                    int *jmax,                 /* number of cells y-direction*/
                    double *alpha,              /* uppwind differencing factor*/
                    double *omg,                /* relaxation factor */
                    double *tau,                /* safety factor for time step*/
                    int *itermax,              /* max. number of iterations  */
                    double *eps,                /* accuracy bound for pressure*/
                    double *dt_value,           /* time for output */
                    char *problem,       /* problem string */
                    char *geometry)    /* path/filename to geometry file */
{
    READ_DOUBLE(szFileName, *xlength);
    READ_DOUBLE(szFileName, *ylength);
    
    READ_DOUBLE(szFileName, *Re);
    READ_DOUBLE(szFileName, *t_end);
    READ_DOUBLE(szFileName, *dt);
    
    READ_INT   (szFileName, *imax);
    READ_INT   (szFileName, *jmax);
    
    READ_DOUBLE(szFileName, *omg);
    READ_DOUBLE(szFileName, *eps);
    READ_DOUBLE(szFileName, *tau);
    READ_DOUBLE(szFileName, *alpha);
    
    READ_INT   (szFileName, *itermax);
    READ_DOUBLE(szFileName, *dt_value);
    
    READ_DOUBLE(szFileName, *UI);
    READ_DOUBLE(szFileName, *VI);
    READ_DOUBLE(szFileName, *GX);
    READ_DOUBLE(szFileName, *GY);
    READ_DOUBLE(szFileName, *PI);
    
    READ_STRING(szFileName, problem);
    READ_STRING(szFileName, geometry);
    
    *dx = *xlength / (double) (*imax);
    *dy = *ylength / (double) (*jmax);
    
    return 1;
}

void init_uvp(
        double UI,
        double VI,
        double PI,
        int imax,
        int jmax,
        double **U,
        double **V,
        double **P
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
        int **Flag
)
{
    int **pic = NULL;
    
    pic = read_pgm(geometry); // NOTE: this is covering just the inner part of the image, so it is imax*jmax
    
    // Set the corner cases
    // bottom-left
    Flag[0][0] = 1;
    // top-left
    Flag[0][jmax + 1] = 1;
    // bottom-right
    Flag[imax + 1][0] = 1;
    // top-right
    Flag[imax + 1][jmax + 1] = 1;
    // Set the outer boundary + the first inner layers
    for (int i = 1; i < imax + 1; ++i)
    {
        // Outer boundary
        Flag[i][0] = 1;
        Flag[i][jmax + 1] = 1;
        // First inner rows
        Flag[i][1] = pic[i - 1][0];
        Flag[i][jmax] = pic[i - 1][jmax - 1];
    }
    for (int j = 1; j < imax + 1; ++j)
    {
        // Outer boundary
        Flag[0][j] = 1;
        Flag[imax + 1][j] = 1;
        // First inner columns
        Flag[1][j] = pic[0][j - 1];
        Flag[imax][j] = pic[imax - 1][j - 1];
    }
    // Set the inner domain geometry
    for (int i = 2; i < imax; i++)
    {
        for (int j = 2; j < jmax; j++)
        {
            Flag[i][j] = pic[i - 1][j - 1];
        }
    }
    // Set the inner domain flags
    for (int i = 1; i < imax + 1; i++)
    {
        for (int j = 1; j < jmax + 1; j++)
        {
            Flag[i][j] += TOP * isObstacle(Flag[i][j + 1])
                          + BOT * isObstacle(Flag[i][j - 1])
                          + LEFT * isObstacle(Flag[i - 1][j])
                          + RIGHT * isObstacle(Flag[i + 1][j]);
            printf("%d ", Flag[i][j]);
        }
        printf("\n");
    }
    
    free_imatrix(pic, 0, imax + 1, 0, jmax + 1);
}
