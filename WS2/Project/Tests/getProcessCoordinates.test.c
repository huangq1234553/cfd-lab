// Testing the uc_comm() function

#include "../parallel.h"
#include "../helper.h"
#include "getProcessCoordinates.test.h"
#include "testing.h"
#include "../logger.h"

static char *TEST_NAME = "getProcessCoordinates.test";

static int processesPerRow, processesPerColumn;
static int omegaIExpected, omegaJExpected; // The coordinates

static void setup(int processesPerRow_, int processesPerColumn_, int omegaI_, int omegaJ_)
{
    processesPerRow = processesPerRow_;
    processesPerColumn = processesPerColumn_;
    omegaIExpected = omegaI_;
    omegaJExpected = omegaJ_;
}

static void teardown()
{
}

int getProcessCoordinatesTest(int mpiRank, int mpiNumProc)
{
    // Now assume we are in a 2x2 grid, so let's check we have enough MPI processes
    assert(mpiNumProc >= 4);
    // Now setup the expectations
    switch (mpiRank)
    {
        case 0:
            setup(2, 2, 0, 0);
            break;
        case 1:
            setup(2, 2, 1, 0);
            break;
        case 2:
            setup(2, 2, 0, 1);
            break;
        case 3:
            setup(2, 2, 1, 1);
            break;
        default:
            return 0; // In case of extra processors, they don't have anything to do
    }
    int omegaI, omegaJ;
    // Perform coordinates computation
    getProcessCoordinates(processesPerRow, mpiRank, &omegaI, &omegaJ);
    
    // Now check values
    int failedTests = 0;
    failedTests += expectEqual(omegaI, omegaIExpected, TEST_NAME, "[R%d] I coordinate check failed: omegaI = ");
    failedTests += expectEqual(omegaJ, omegaJExpected, TEST_NAME, "[R%d] J coordinate check failed: omegaJ = ");
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
