#include <mpi.h>
#include <stdio.h>
#include "../logger.h"
#include "uvComm.test.h"

/// Setup and teardown functions
void setupGlobal(int *argc, char ***argv, int *mpiRank, int *mpiNumProc)
{
    // Setup MPI
    MPI_Status status;
    MPI_Init(argc, argv);
    MPI_Comm_size(MPI_COMM_WORLD, mpiNumProc);
    MPI_Comm_rank(MPI_COMM_WORLD, mpiRank);
    // Setup logging
    setLoggerDebugLevel(INFO);
    char fnamebuf[128];
    sprintf(fnamebuf, "tests.%d.log", *mpiRank);
    setLoggerFileName(fnamebuf); // It will log in current working directory
    setLoggerStartTime();
    openLogFile();
}

void teardownGlobal()
{
    // Teardown logging
    closeLogFile();
    // Teardown MPI
    MPI_Finalize();
}

// Main
int main(int argc, char **argv)
{
    int mpiRank, mpiNumProc;
    setupGlobal(&argc, &argv, &mpiRank, &mpiNumProc);
    // Exec tests
    if (mpiRank == 0)
    {
        logMsg("Starting tests! :)");
    }
    int uvCommFailures = uvCommTest(mpiRank, mpiNumProc);

//  Now gather the sum of all failed tests
    int globalUvCommFailures = 0;
    MPI_Reduce(&uvCommFailures, &globalUvCommFailures, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    teardownGlobal();
    return (globalUvCommFailures >= 0);
}

//eof
