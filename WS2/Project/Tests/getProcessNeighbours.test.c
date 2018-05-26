// Testing the uc_comm() function

#include "../parallel.h"
#include "../helper.h"
#include "getProcessNeighbours.test.h"
#include "testing.h"
#include "../logger.h"

static char *TEST_NAME = "getProcessNeighbours.test";

static int processesPerRow, processesPerColumn, omegaI, omegaJ;
static int rankLExpected, rankRExpected, rankBExpected, rankTExpected;

static void
setup(int processesPerRow_, int processesPerColumn_, int omegaI_, int omegaJ_, int rank_l_, int rank_r_, int rank_b_,
      int rank_t_)
{
    processesPerRow = processesPerRow_;
    processesPerColumn = processesPerColumn_;
    omegaI = omegaI_;
    omegaJ = omegaJ_;
    rankLExpected = rank_l_;
    rankRExpected = rank_r_;
    rankBExpected = rank_b_;
    rankTExpected = rank_t_;
}

static void teardown()
{
}

int getProcessNeighboursTest(int mpiRank, int mpiNumProc)
{
    // Now assume we are in a 2x2 grid, so let's check we have enough MPI processes
    assert(mpiNumProc >= 4);
    // Now setup the expectations
    switch (mpiRank)
    {
        case 0:
            setup(2, 2, 0, 0, MPI_PROC_NULL, 1, MPI_PROC_NULL, 2);
            break;
        case 1:
            setup(2, 2, 1, 0, 0, MPI_PROC_NULL, MPI_PROC_NULL, 3);
            break;
        case 2:
            setup(2, 2, 0, 1, MPI_PROC_NULL, 3, 0, MPI_PROC_NULL);
            break;
        case 3:
            setup(2, 2, 1, 1, 2, MPI_PROC_NULL, 1, MPI_PROC_NULL);
            break;
        default:
            return 0; // In case of extra processors, they don't have anything to do
    }
    int rankL, rankR, rankB, rankT;
    // Perform coordinates computation
    getProcessNeighbours(processesPerRow, processesPerColumn, &omegaI, &omegaJ, &rankL, &rankR, &rankB, &rankT);
    
    // Now check values
    int failedTests = 0;
    failedTests += expectEqual(rankL, rankLExpected, TEST_NAME, "[R%d] Left neighbor check failed: rankL = ");
    failedTests += expectEqual(rankR, rankRExpected, TEST_NAME, "[R%d] Right neighbor check failed: rankR = ");
    failedTests += expectEqual(rankB, rankBExpected, TEST_NAME, "[R%d] Bottom neighbor check failed: rankB = ");
    failedTests += expectEqual(rankT, rankTExpected, TEST_NAME, "[R%d] Top neighbor check failed: rankT = ");
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
