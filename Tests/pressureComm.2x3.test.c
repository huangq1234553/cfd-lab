#include "pressureComm.test.h"
#include "../parallel.h"
#include "../logger.h"
#include "../helper.h"
#include "testing.h"

static char *TEST_NAME = "pressureComm.2x3.test";

static double **P, *bufSend, *bufRecv;
static int rank_l, rank_r, rank_b, rank_t;
static int imax = 2, jmax = 2;
static MPI_Status *status;

static void setup(int mpiRank, int rank_l_, int rank_r_, int rank_b_, int rank_t_)
{
    logTestEvent(DEBUG, TEST_NAME, "Setup start...");
    P = matrix(0, imax + 2, 0, jmax + 2);
    init_matrix(P, 0, imax + 2, 0, jmax + 2, mpiRank); // Initialize with rank
    bufSend = (double *) malloc((size_t) (3 * (imax + 3) * sizeof(double *)));
    bufRecv = (double *) malloc((size_t) (3 * (imax + 3) * sizeof(double *)));
    // rank
    rank_l = rank_l_;
    rank_r = rank_r_;
    rank_b = rank_b_;
    rank_t = rank_t_;
    logTestEvent(DEBUG, TEST_NAME, "Setup end!");
}
static void teardown()
{
    free_matrix(P, 0, imax + 2, 0, jmax + 2);
    free(bufSend);
    free(bufRecv);
}

//
int pressureCommTest23(int mpiRank, int mpiNumProc)
{
    // Now assume we are in a 2x2 grid, so let's check we have enough MPI processes
    assert(mpiNumProc >= 6);
    // Barrier to make sure all processes are in sync here
    MPI_Barrier(MPI_COMM_WORLD);
    // Now setup neighbors
    switch (mpiRank)
    {
        case 0:
            setup(mpiRank, MPI_PROC_NULL, 1, MPI_PROC_NULL, 2);
            break;
        case 1:
            setup(mpiRank, 0, MPI_PROC_NULL, MPI_PROC_NULL, 3);
            break;
        case 2:
            setup(mpiRank, MPI_PROC_NULL, 3, 0, 4);
            break;
        case 3:
            setup(mpiRank, 2, MPI_PROC_NULL, 1, 5);
            break;
        case 4:
            setup(mpiRank, MPI_PROC_NULL, 5, 2, MPI_PROC_NULL);
            break;
        case 5:
            setup(mpiRank, 4, MPI_PROC_NULL, 3, MPI_PROC_NULL);
            break;
        default:
            MPI_Barrier(MPI_COMM_WORLD); // Sync at exit
            return 0; // In case of extra processors, they don't have anything to do
    }
    // Perform communication
    pressure_comm(P, rank_l, rank_r, rank_b, rank_t, bufSend, bufRecv, status, imax, jmax);
    
    // Now check values
    int failedTests = 0;
    if (rank_l != MPI_PROC_NULL) //left
    {
        for (int i = 1; i <= imax; ++i)
        {
            failedTests += expectEqual(P[0][i], rank_l, TEST_NAME, "[R%d] Left ghost layer check failed: P[0][%d] = ",
                                       mpiRank, i);
        }
    }
    
    if (rank_r != MPI_PROC_NULL) //right
    {
        for (int i = 1; i <= imax; ++i)
        {
            logTestEvent(DEBUG, TEST_NAME, "[R%d] Checking right ghost layer at %d", mpiRank, i);
            failedTests += expectEqual(P[imax + 1][i], rank_r, TEST_NAME,
                                       "[R%d] Right ghost layer check failed: P[imax+1][%d] = ", mpiRank, i);
        }
    }
    
    if (rank_b != MPI_PROC_NULL) //bottom
    {
        for (int i = 1; i <= imax; ++i)
        {
            failedTests += expectEqual(P[i][0], rank_b, TEST_NAME, "[R%d] Bottom ghost layer check failed: P[%d][0] = ",
                                       mpiRank, i);
        }
    }
    
    if (rank_t != MPI_PROC_NULL) //top
    {
        for (int i = 1; i <= imax; ++i)
        {
            failedTests += expectEqual(P[i][jmax + 1], rank_t, TEST_NAME,
                                       "[R%d] Top ghost layer check failed: P[%d][jmax+1] = ", mpiRank, i);
        }
    }
    //
    
    if (failedTests == 0)
    {
        logTestEvent(INFO, TEST_NAME, "[R%d] Test completed successfully! :)", mpiRank);
    }
    else
    {
        logTestEvent(INFO, TEST_NAME, "[R%d] Some tests failed (%d) :(", mpiRank, failedTests);
    }
    // Teardown
    teardown();
    
    // Barrier to make sure all processes are in sync here
    MPI_Barrier(MPI_COMM_WORLD);
    
    return failedTests;
}

// eof
