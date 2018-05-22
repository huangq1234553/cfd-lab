#include "parallel.h"


void Program_Message(char *txt)
/* produces a stderr text output  */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
}


void Programm_Sync(char *txt)
/* produces a stderr textoutput and synchronize all processes */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                             /* synchronize output */  
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
}


void Programm_Stop(char *txt)
/* all processes will produce a text output, be synchronized and finished */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                           /* synchronize output */
   fprintf(stderr,"-STOP- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   exit(1);
}


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
   int imax,
   int jmax){

   // Send right || Receive left

   // Fill bufSend with second-rightmost column for sending
   for(int j = 1; j < jmax+1; ++j){
      bufSend[j-1] = P[imax][j];
   }
   
   MPI_Sendrecv(
    bufSend,                                      // send buffer
    jmax,               // number of records
    MPI_FLOAT,                                   // datatype
    // destination process
    rank_r,
    1,                                           // message tag
    bufRecv,                                  // receive buffer
    jmax,               // number of records
    MPI_FLOAT,                                   // datatype
    // source process
    rank_l,
    1,                                           // message tag
    MPI_COMM_WORLD,                              // communicator
    status                                      // status of communication
  );

   // Fill leftmost column in P for bufRecv
   for(int j = 1; j < jmax+1; ++j){
      P[0][j] = bufRecv[j-1];
   }
   
   // Send left || Receive right
   
   // Fill bufSend with second-leftmost column for sending
   for(int j = 1; j < jmax+1; ++j){
      bufSend[j-1] = P[1][j];
   }

   MPI_Sendrecv(
    bufSend,                                      // send buffer
    jmax,               // number of records
    MPI_FLOAT,                                   // datatype
    rank_l,                                         // destination process
    2,                                           // message tag
    bufRecv,                                  // receive buffer
    jmax,               // number of records
    MPI_FLOAT,                                   // datatype
    rank_r,                                 // source process
    2,                                           // message tag
    MPI_COMM_WORLD,                              // communicator
    status                                      // status of communication
  );

   // Fill rightmost column of P with bufRecv
   for(int j = 1; j < jmax + 1; ++j){
      P[imax +1][j] = bufRecv[j-1];
   }


   // Send bottom || Receive top   

   // Fill bufSend with second-bottommost row for sending
   for(int i = 1; i < imax + 1; ++i){
      bufSend[i-1] = P[i][1];
   }

   MPI_Sendrecv(
    bufSend,                                      // send buffer
    imax,               // number of records
    MPI_FLOAT,                                   // datatype
    rank_b,                                 // destination process
    3,                                           // message tag
    bufRecv,                                  // receive buffer
    imax,               // number of records
    MPI_FLOAT,                                   // datatype
    rank_t,                                     // source process
    3,                                           // message tag
    MPI_COMM_WORLD,                              // communicator
    status                                      // status of communication
  );

   // Fill topmost row in P with bufRecv
   for(int i = 1; i < imax + 1; ++i){
      P[i][jmax+1] = bufRecv[i-1];
   }

   // Send top || Receive bottom

   // Fill bufSend with second-topmost column for sending
   for(int i = 1; i < imax + 1; ++i){
      bufSend[i-1] = P[i][jmax];
   }

   MPI_Sendrecv(
    bufSend,                                      // send buffer
    (ir - il) - 1,               // number of records
    MPI_FLOAT,                                   // datatype
    rank_t,                                 // destination process
    4,                                           // message tag
    bufRecv,                                  // receive buffer
    (ir - il) - 1,               // number of records
    MPI_FLOAT,                                   // datatype
    rank_b,                                      // source process
    4,                                           // message tag
    MPI_COMM_WORLD,                              // communicator
    status                                      // status of communication
  );

   // Fill second-topmost column in P with bufRecv
   for(int i = 1; i < imax + 1; ++i){
      P[i][0] = bufRecv[i-1];
   }

}

void init_parallel (
int iproc,int jproc,int imax,int jmax,
int myrank,int *il,int *ir,int *jb,int *jt,
int *rank_l,int *rank_r,int *rank_b,int *rank_t,
int *omg_i,int *omg_j,int num_proc){

  int imax_local, jmax_local;
  // Starting (omg_i,omg_j) indices form (0,0)

  if(imax%iproc != 0 && jmax%jproc != 0){
    printf("This geometry is not allowed, please choose a suitable number of processes or change grid size");
  }
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

  (*il) = imax_local*(*omg_i);
  (*ir) = imax_local*(*omg_i+1) - 1;
  (*jb) = jmax_local*(*omg_j);
  (*jt) = jmax_local*(*omg_j+1) - 1;

}
