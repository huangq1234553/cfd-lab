#include "parallel.h"


void Program_Message(char *txt)
/* produces a stderr text output  */

{
    int myrank;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    fprintf(stderr, "-MESSAGE- P:%2d : %s\n", myrank, txt);
    fflush(stdout);
    fflush(stderr);
}


void Programm_Sync(char *txt)
/* produces a stderr textoutput and synchronize all processes */

{
    int myrank;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Barrier(MPI_COMM_WORLD);                             /* synchronize output */
    fprintf(stderr, "-MESSAGE- P:%2d : %s\n", myrank, txt);
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
    fprintf(stderr, "-STOP- P:%2d : %s\n", myrank, txt);
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
        int jmax)
{
    
    // Send right || Receive left
    
    // Fill bufSend with second-rightmost column for sending
    for (int j = 1; j < jmax + 1; ++j)
    {
        bufSend[j - 1] = P[imax][j];
    }
    
    MPI_Sendrecv(
            bufSend,                                      // send buffer
            jmax,                           // number of records
            MPI_FLOAT,                                   // datatype
            // destination process
            rank_r,
            1,                                           // message tag
            bufRecv,                                  // receive buffer
            jmax,                               // number of records
            MPI_FLOAT,                                   // datatype
            // source process
            rank_l,
            1,                                           // message tag
            MPI_COMM_WORLD,                              // communicator
            status                                      // status of communication
    );


    
    // Fill leftmost column in P for bufRecv
    
    if (rank_l != MPI_PROC_NULL)
    {
        for (int j = 1; j < jmax + 1; ++j)
        {
            P[0][j] = bufRecv[j - 1];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // Send left || Receive right
    
    // Fill bufSend with second-leftmost column for sending
    for (int j = 1; j < jmax + 1; ++j)
    {
        bufSend[j - 1] = P[1][j];
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
    if (rank_r != MPI_PROC_NULL)
    {
        for (int j = 1; j < jmax + 1; ++j)
        {
            P[imax + 1][j] = bufRecv[j - 1];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Send bottom || Receive top
    
    // Fill bufSend with second-bottommost row for sending
    for (int i = 1; i < imax + 1; ++i)
    {
        bufSend[i - 1] = P[i][1];
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
    if (rank_t != MPI_PROC_NULL)
    {
        for (int i = 1; i < imax + 1; ++i)
        {
            P[i][jmax + 1] = bufRecv[i - 1];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // Send top || Receive bottom
    
    // Fill bufSend with second-topmost column for sending
    for (int i = 1; i < imax + 1; ++i)
    {
        bufSend[i - 1] = P[i][jmax];
    }
    
    MPI_Sendrecv(
            bufSend,                                      // send buffer
            imax,               // number of records
            MPI_FLOAT,                                   // datatype
            rank_t,                                 // destination process
            4,                                           // message tag
            bufRecv,                                  // receive buffer
            imax,               // number of records
            MPI_FLOAT,                                   // datatype
            rank_b,                                      // source process
            4,                                           // message tag
            MPI_COMM_WORLD,                              // communicator
            status                                      // status of communication
    );
    
    // Fill second-topmost column in P with bufRecv
    if (rank_b != MPI_PROC_NULL)
    {
        for (int i = 1; i < imax + 1; ++i)
        {
            P[i][0] = bufRecv[i - 1];
        }
    }
    
}


void uv_comm(double **U, double **V, int rank_l, int rank_r, int rank_b, int rank_t,
             double *bufSend, double *bufRecv, MPI_Status *status, int imax, int jmax)
{
    
    // Send right || Receive left
    int count = 0;
    // Fill bufSend with third-rightmost column of U
    for (int j = 1; j < jmax + 3; ++j)
    {
        bufSend[j - 1] = U[imax + 3 - 2][j];
        count++;
    }
    // Fill bufSend with second-rightmost column of V
    for (int j = 1; j < jmax + 3; ++j)
    {
        bufSend[(j - 1) + (jmax + 2)] = V[imax + 2 - 1][j];
        count++;
    }
    
    
    MPI_Sendrecv(
            bufSend,                                      // send buffer
            count,               // number of records
            MPI_DOUBLE,                                   // datatype
            // destination process
            rank_r,
            1,                                           // message tag
            bufRecv,                                  // receive buffer
            count,               // number of records
            MPI_DOUBLE,                                   // datatype
            // source process
            rank_l,
            1,                                           // message tag
            MPI_COMM_WORLD,                              // communicator
            status                                      // status of communication
    );
    
    if (rank_l != MPI_PROC_NULL)
    {
        // Fill first-leftmost column of U
        for (int j = 1; j < jmax + 3; ++j)
        {
            U[0][j] = bufRecv[j - 1];
        }
        // Fill first-leftmost column of V
        for (int j = 1; j < jmax + 3; ++j)
        {
            V[0][j] = bufRecv[(j - 1) + (jmax + 2)];
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    // Send left || Receive right
    
    count = 0;
    // Fill bufSend with third-leftmost column of U
    for (int j = 1; j < jmax + 2; ++j)
    {
        bufSend[j - 1] = U[2][j];
        count++;
    }
    // Fill bufSend with second-leftmost column of V
    for (int j = 1; j < jmax + 3; ++j)
    {
        bufSend[(j - 1) + (jmax + 1)] = V[1][j];
        count++;
    }
    
    MPI_Sendrecv(
            bufSend,                                      // send buffer
            count,               // number of records
            MPI_FLOAT,                                   // datatype
            rank_l,                                         // destination process
            2,                                           // message tag
            bufRecv,                                  // receive buffer
            count,               // number of records
            MPI_FLOAT,                                   // datatype
            rank_r,                                 // source process
            2,                                           // message tag
            MPI_COMM_WORLD,                              // communicator
            status                                      // status of communication
    );
    
    if (rank_r != MPI_PROC_NULL)
    {
        // Fill first-rightmost column of U
        for (int j = 1; j < jmax + 2; ++j)
        {
            U[imax + 3][j] = bufRecv[j - 1];
        }
        // Fill first-rightmost column of V
        for (int j = 1; j < jmax + 3; ++j)
        {
            V[imax + 2][j] = bufRecv[(j - 1) + (jmax + 1)];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // Send bottom || Receive top
    
    count = 0;
    // Fill bufSend with second-bottommost row of U
    for (int i = 1; i < imax + 3; ++i)
    {
        bufSend[i - 1] = U[i][1];
        count++;
    }
    // Fill bufSend with third-bottommost row of V
    for (int i = 1; i < imax + 2; ++i)
    {
        bufSend[(i - 1) + (imax + 2)] = V[i][2];
        count++;
    }
    
    MPI_Sendrecv(
            bufSend,                                      // send buffer
            count,               // number of records
            MPI_FLOAT,                                   // datatype
            rank_b,                                 // destination process
            3,                                           // message tag
            bufRecv,                                  // receive buffer
            count,               // number of records
            MPI_FLOAT,                                   // datatype
            rank_t,                                     // source process
            3,                                           // message tag
            MPI_COMM_WORLD,                              // communicator
            status                                      // status of communication
    );
    
    if (rank_t != MPI_PROC_NULL)
    {
        // Fill first-topmost row of U
        for (int i = 1; i < imax + 3; ++i)
        {
            U[i][jmax + 2] = bufRecv[i - 1];
        }
        // Fill first-topmost row of V
        for (int i = 1; i < imax + 2; ++i)
        {
            V[i][imax + 3] = bufRecv[(i - 1) + (imax + 2)];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // Send top || Receive bottom
    
    count = 0;
    // Fill bufSend with second-topmost row of U
    for (int i = 1; i < imax + 3; ++i)
    {
        bufSend[i - 1] = U[i][jmax + 2 - 1];
        count++;
    }
    // Fill bufSend with third-topmost row of V
    for (int i = 1; i < imax + 3; ++i)
    {
        bufSend[(i - 1) + (imax + 2)] = V[i][jmax + 3 - 2];
        count++;
    }
    
    MPI_Sendrecv(
            bufSend,                                      // send buffer
            count,               // number of records
            MPI_FLOAT,                                   // datatype
            rank_t,                                 // destination process
            4,                                           // message tag
            bufRecv,                                  // receive buffer
            count,               // number of records
            MPI_FLOAT,                                   // datatype
            rank_b,                                      // source process
            4,                                           // message tag
            MPI_COMM_WORLD,                              // communicator
            status                                      // status of communication
    );
    
    if (rank_b != MPI_PROC_NULL)
    {
        // Fill first-bottommost row of U
        for (int i = 1; i < imax + 3; ++i)
        {
            U[i][0] = bufRecv[i - 1];
        }
        // Fill first-bottommost row of V
        for (int i = 1; i < imax + 3; ++i)
        {
            V[i][0] = bufRecv[(i - 1) + (imax + 2)];
        }
    }
    
}

void init_parallel(
        int iproc, int jproc, int imax, int jmax,
        int myrank, int *il, int *ir, int *jb, int *jt,
        int *rank_l, int *rank_r, int *rank_b, int *rank_t,
        int *omg_i, int *omg_j, int num_proc)
{
    int imax_local, jmax_local;
    // Starting (omg_i,omg_j) indices form (0,0)
    
    if (imax % iproc != 0 && jmax % jproc != 0)
    {
        printf("This geometry is not allowed, please choose a suitable number of processes or change grid size");
    }
    *omg_i = myrank % iproc;
    *omg_j = myrank / iproc;
    
    // Here we try to determine the ranks of our neighbours (l/r/b/t)
    if ((*omg_i != 0))
    {
        (*rank_l) = (*omg_i - 1) + (*omg_j) * iproc;
    }
    else
    {
        (*rank_l) = MPI_PROC_NULL;
    }
    
    
    if ((*omg_i != iproc - 1))
    {
        (*rank_r) = (*omg_i + 1) + (*omg_j) * iproc;
    }
    else
    {
        (*rank_r) = MPI_PROC_NULL;
    }
    
    if ((*omg_j != 0))
    {
        (*rank_b) = (*omg_i) + (*omg_j - 1) * iproc;
    }
    else
    {
        (*rank_b) = MPI_PROC_NULL;
    }
    
    
    if ((*omg_j != jproc - 1))
    {
        (*rank_t) = (*omg_i) + (*omg_j + 1) * iproc;
    }
    else
    {
        (*rank_t) = MPI_PROC_NULL;
    }
    
    // Assume that imax%iproc == 0 && jmax%jproc == 0
    imax_local = imax / iproc;
    jmax_local = jmax / jproc;
    
    (*il) = imax_local * (*omg_i);
    (*ir) = imax_local * (*omg_i + 1) - 1;
    (*jb) = jmax_local * (*omg_j);
    (*jt) = jmax_local * (*omg_j + 1) - 1;
    
}
