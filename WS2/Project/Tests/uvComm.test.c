// Testing the uc_comm() function

#include "../parallel.h"
#include "../helper.h"
#include "uvComm.test.h"
#include "testing.h"
#include "../logger.h"

static char *TEST_NAME = "uvComm.test";

static double **U, **V, *bufSend, *bufRecv;
static int rank_l, rank_r, rank_b, rank_t;
static int imax = 2, jmax = 2;
static MPI_Status *status;

// Here we need to declare as static to make the function local in scope
static void setup(int mpiRank, int rank_l_, int rank_r_, int rank_b_, int rank_t_)
{
    logTestEvent(DEBUG, TEST_NAME, "Setup start...");
    U = matrix(0, imax + 3, 0, jmax + 2);
    V = matrix(0, imax + 2, 0, jmax + 3);
    init_matrix(U, 0, imax + 3, 0, jmax + 2, mpiRank); // Initialize with rank
    init_matrix(V, 0, imax + 2, 0, jmax + 3, mpiRank); // Initialize with rank
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
    free_matrix(U, 0, imax + 2, 0, jmax + 1);
    free_matrix(V, 0, imax + 1, 0, jmax + 2);
    free(bufSend);
    free(bufRecv);
}

int uvCommTest(int mpiRank, int mpiNumProc)
{
    // Now assume we are in a 2x2 grid, so let's check we have enough MPI processes
    assert(mpiNumProc >= 4);
//    assertEqual(mpiNumProc,4);
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
            setup(mpiRank, MPI_PROC_NULL, 3, 0, MPI_PROC_NULL);
            break;
        case 3:
            setup(mpiRank, 2, MPI_PROC_NULL, 1, MPI_PROC_NULL);
            break;
        default:
            return 0; // In case of extra processors, they don't have anything to do
    }
    // Perform communication
    uv_comm(U, V, rank_l, rank_r, rank_b, rank_t, bufSend, bufRecv, status, imax, jmax);
    
    // Now check values
    int failedTests = 0;
    if (rank_l != MPI_PROC_NULL) //left
    {
        for (int i = 1; i <= imax; ++i)
        {
            failedTests += expectEqual(U[0][i], rank_l, TEST_NAME, "[R%d] Left ghost layer check failed: U[0][%d] = ",
                                       mpiRank, i);
            failedTests += expectEqual(V[0][i], rank_l, TEST_NAME, "[R%d] Left ghost layer check failed: V[0][%d] = ",
                                       mpiRank, i);
        }
    }
    
    if (rank_r != MPI_PROC_NULL) //right
    {
        for (int i = 1; i <= imax; ++i)
        {
            logTestEvent(DEBUG, TEST_NAME, "[R%d] Checking right ghost layer at %d", mpiRank, i);
            failedTests += expectEqual(U[imax + 3][i], rank_r, TEST_NAME,
                                       "[R%d] Right ghost layer check failed: U[imax+3][%d] = ", mpiRank, i, U);
            failedTests += expectEqual(V[imax + 2][i], rank_r, TEST_NAME,
                                       "[R%d] Right ghost layer check failed: V[imax+2][%d] = ", mpiRank, i, V);
        }
    }
    
    if (rank_b != MPI_PROC_NULL) //bottom
    {
        for (int i = 1; i <= imax; ++i)
        {
            failedTests += expectEqual(U[i][0], rank_b, TEST_NAME, "[R%d] Bottom ghost layer check failed: U[%d][0] = ",
                                       mpiRank, i);
            failedTests += expectEqual(V[i][0], rank_b, TEST_NAME, "[R%d] Bottom ghost layer check failed: V[%d][0] = ",
                                       mpiRank, i);
        }
    }
    
    if (rank_t != MPI_PROC_NULL) //top
    {
        for (int i = 1; i <= imax; ++i)
        {
            failedTests += expectEqual(U[i][jmax + 2], rank_t, TEST_NAME,
                                       "[R%d] Top ghost layer check failed: U[%d][jmax+2] = ", mpiRank, i, U);
            failedTests += expectEqual(V[i][jmax + 3], rank_t, TEST_NAME,
                                       "[R%d] Top ghost layer check failed: V[%d][jmax+3] = ", mpiRank, i, V);
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
    
    return failedTests;
}

//eof
