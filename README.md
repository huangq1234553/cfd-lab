# cfd-lab

Repository for the code of GroupA, CFD Lab course SoSe2018 at TUM.
Enjoy! :)

## README - How to execute the code and the test cases

1) Easy way to get started: execute the unit tests and pre-defined simulation case!

    Compile:  
    `$ make`  
    Run unit tests:  
    `$ make runtests`  
    Run simulation with 4 MPI processes:  
    `$ make runsim`  
    Run simulation with sequential code (1 MPI process):  
    `$ make runsimseq`  
    Run simulation with 16 MPI processes:  
    `$ make runsimpar16`  
    The above simulations will create visualization files inside the `Out/` subfolder.
    

2) Manually running (example):

    `$ make`  
    `$ mpirun -np <numProcesses> ./sim -q problem.dat`

## Supported command line arguments (OPTIONS)
List of supported arguments:

    --iproc
        Overrides iproc value in .dat file.

    --jproc
        Overrides jproc value in .dat file.

    -q
        "Quiet" mode: reduces the verbosity of output/log by disabling INFO level traces
        (it still retains PRODUCTION, WARNING, ERROR level traces).

    -o path/to/out/folder/
        Output directory path. In case this is not set, an "Out" subfolder will be created inside the folder containing
        the configuration .dat file.

## Logging
The same output which is written to stdout is also written into a log file (sim.log) in the output folder, along with
the output visualization files.
This is useful to keep trace of what parameters where used to generate the output and to investigate potential issues.

