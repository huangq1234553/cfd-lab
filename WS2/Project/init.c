#include "helper.h"
#include "init.h"

int read_parameters(const char *szFileName, double *Re, double *UI, double *VI, double *PI, double *GX, double *GY,
                    double *t_end, double *xlength, double *ylength, double *dt, double *dx, double *dy, int *imax,
                    int *jmax, double *alpha, double *omg, double *tau, int *itermax, double *eps, double *dt_value,
                    char *problem, char *geometry, int *iproc, int *jproc)    /* path/filename to geometry file */
{
    READ_DOUBLE( szFileName, *xlength );
    READ_DOUBLE( szFileName, *ylength );
    
    READ_DOUBLE( szFileName, *Re    );
    READ_DOUBLE( szFileName, *t_end );
    READ_DOUBLE( szFileName, *dt    );
    
    READ_INT   ( szFileName, *imax );
    READ_INT   ( szFileName, *jmax );

    READ_INT   ( szFileName, *iproc );
    READ_INT   ( szFileName, *jproc );
    
    READ_DOUBLE( szFileName, *omg   );
    READ_DOUBLE( szFileName, *eps   );
    READ_DOUBLE( szFileName, *tau   );
    READ_DOUBLE( szFileName, *alpha );
    
    READ_INT   ( szFileName, *itermax );
    READ_DOUBLE( szFileName, *dt_value );
    
    READ_DOUBLE( szFileName, *UI );
    READ_DOUBLE( szFileName, *VI );
    READ_DOUBLE( szFileName, *GX );
    READ_DOUBLE( szFileName, *GY );
    READ_DOUBLE( szFileName, *PI );
    
    // READ_STRING( szFileName, problem);
    // READ_STRING( szFileName, geometry);

    READ_INT   ( szFileName, *iproc );
    READ_INT   ( szFileName, *jproc );

   *dx = *xlength / (double)(*imax);
   *dy = *ylength / (double)(*jmax);

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
){

  init_matrix( P, 0, imax + 1   , 0   , jmax + 1, PI);
  init_matrix( U, 0, imax + 2   , 0   , jmax + 1, UI);
  init_matrix( V, 0, imax + 1   , 0   , jmax + 2, VI);
  
}




