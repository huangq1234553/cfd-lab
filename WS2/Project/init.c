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
    
    READ_STRING( szFileName, problem);
    READ_STRING( szFileName, geometry);

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
  init_matrix( U, 0, imax+1, 0, jmax+1, UI);
  init_matrix( V, 0, imax+1, 0, jmax+1, VI);
  init_matrix( P, 0, imax+1, 0, jmax+1, PI);
}



void init_parallel (
int iproc,int jproc,int imax,int jmax,
int myrank,int *il,int *ir,int *jb,int *jt,
int *rank_l,int *rank_r,int *rank_b,int *rank_t,
int *omg_i,int *omg_j,int num_proc){

  int imax_local, jmax_local;
  // Starting (omg_i,omg_j) indices form (0,0)

  *omg_i = myrank%iproc;
  *omg_j = myrank/iproc; 

  if((*omg_i != 0)){
    (*rank_l) = (*omg_i-1) + (*omg_j)*iproc;
  }
  else{
    (*rank_l) = MPI_PROC_NULL;
  }


  if((*omg_i != iproc-1)){
    (*rank_r) = (*omg_i+1) + (*omg_j)*iproc;
  }
  else{
    (*rank_r) = MPI_PROC_NULL;
  }
  
  if((*omg_j != 0)){
    (*rank_b) = (*omg_i) + (*omg_j-1)*iproc;
  }
  else{
    (*rank_b) = MPI_PROC_NULL;
  }


  if((*omg_j != jproc-1)){
    (*rank_t) = (*omg_i) + (*omg_j+1)*iproc;
  }
  else{
    (*rank_t) = MPI_PROC_NULL;
  }

  // Assume that imax%iproc == 0 && jmax%jproc == 0
  imax_local = imax/iproc;
  jmax_local = jmax/jproc;

  (*il) = imax_local*(omg_i);
  (*ir) = imax_local*(omg_i+1) - 1;
  (*jb) = jmax_local*(omg_j);
  (*jt) = jmax_local*(omg_j+1) - 1;

}

