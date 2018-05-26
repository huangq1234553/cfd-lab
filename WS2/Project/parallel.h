#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

void Program_Message(char *txt);
/* produces a stderr text output  */



void Programm_Sync(char *txt);
/* produces a stderr textoutput and synchronize all processes */



void Programm_Stop(char *txt);
/* all processes will produce a text output, be synchronized and finished */

void pressure_comm(
   double **P,
   int il,
   int ir,
   int jb,
   int jt,
   int rank_l,
   int rank_r,
   int rank_b,
   int rank_t,
   double *bufSend,
   double *bufRecv,
   MPI_Status *status,
   int imax_local,
   int jmax_local);

void uv_comm( double**U,double**V, int rank_l,int rank_r,int rank_b,int rank_t,
              double *bufSend,double *bufRecv, MPI_Status *status, int imax, int jmax);

void init_parallel (
  int iproc,
  int jproc,
  int imax,
  int jmax,
  int myrank,
  int *il,
  int *ir,
  int *jb,
  int *jt,
  int *rank_l,
  int *rank_r,
  int *rank_b,
  int *rank_t,
  int *omg_i,
  int *omg_j,
  int num_proc);

void getProcessCoordinates(int iproc, int myrank, int *omg_i, int *omg_j);

void getProcessNeighbours(int iproc, int jproc, const int *omg_i, const int *omg_j, int *rank_l, int *rank_r,
                          int *rank_b,
                          int *rank_t);
