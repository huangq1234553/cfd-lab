# cfd-lab

Repository for the code of GroupA, CFD Lab course SoSe2018 at TUM.
Enjoy! :)

Please download the attached videos for the various cases from our github if you don't want to run for hours:
https://github.com/huangq1234553/cfd-lab/tree/project

## README - How to execute the code and the test cases

1) Easy way to get started: execute the pre-defined test cases!

    `$ make`  
2) Run the code for the various cases:

	Basic flow over step. Modified to run with extreme element removal

    $ ./sim ./Cases/FlowOverStep/FlowOverStep.dat -q
	
	Flow over step with shifted outflow region. Modified to run with low velocities

    $ ./sim ./Cases/ComplexFlowOverStep/ComplexFlowOverStep.dat -q

    Flow over step with shifted obstacle. Additional higher resolution image created.

    $ ./sim ./Cases/Obstacle/ShiftedFlowOverStep.dat -q
    $ ./sim ./Cases/Obstacle/ShiftedFlowOverStep_x4.dat -q

	Obstacle fully immersed in fluid. Additional higher resolution image created.

    $ ./sim ./Cases/Obstacle/Obstacle.dat -q
    $ ./sim ./Cases/Obstacle/Obstacle_x2.dat -q

    Simple bend.
    
    $ ./sim ./Cases/SimpleBend/SimpleBend.dat -q

## Supported command line arguments (Usage & options)
List of supported arguments:

USAGE:

    ./sim configuration.dat [OPTIONS]

OPTIONS:

    -o PATH   Set the output folder path (default is INPUT_DIR/Out).

    --q       Set the logging level to PRODUCTION (ERROR, WARNING and PRODUCTION traces will be enabled).

    --debug   Set the logging level to DEBUG (all traces will be enabled).
  
    --notemp  Disable temperature computation.
  
    --fix-initial-geometry
            Allows auto-fixing the initial geometry before starting the simulation.

INPUT_DIR is the directory where the configuration.dat file is placed.
Default logging level is INFO.
Logging levels ordered by decreasing priority are ERROR, WARNING, PRODUCTION, INFO, DEBUG.

## Logging
The same output which is written to stdout is also written into a log file (sim.log) in the output folder, along with
the output visualization files.
This is useful to keep trace of what parameters where used to generate the output and to investigate potential issues.

## Supported PGM details
Our code supports PGM files in the format as per example published on Moodle, where:
  - the domain outer halo is included in the PGM file
  - 0=NOSLIP, 1=FREESLIP, 2=OUTFLOW, 3=INFLOW, 4=COUPLING, 6=FLUID 

## Configuring parameters in the .dat configuration file
If a variable is omitted in the configuration file, it will take its default value.

  - left_boundary_Dirichlet_U, right_boundary_Dirichlet_U, top_boundary_Dirichlet_U, bottom_boundary_Dirichlet_U  
    left_boundary_Dirichlet_V, right_boundary_Dirichlet_V, top_boundary_Dirichlet_V, bottom_boundary_Dirichlet_V  
    Accepted values: any double  
    Default value: 0.0

  - left_boundary_temp_type, right_boundary_temp_type, top_boundary_temp_type, bottom_boundary_temp_type  
    Accepted values: DIRICHLET, NEUMANN  
    Default value: NEUMANN

  - left_boundary_T, right_boundary_T, top_boundary_T, bottom_boundary_T  
    left_boundary_qN, right_boundary_qN, top_boundary_qN, bottom_boundary_qN  
    Accepted values: any double  
    Default value: 0.0

  - left_boundary_k, right_boundary_k, top_boundary_k, bottom_boundary_k  
    Accepted values: any double  
    Default value: 1.0

For geometry optimization, please use pressure boundary conditions for inflow and outflow regions. To enforce a pressure boundary condition, fill the inflow region with FLUID cells and enforce the following boundary conditions:

   - left_boundary_P, right_boundary_P, top_boundary_P, bottom_boundary_P
     Accepted values: any double
     Default value: 0.0

Other geometry optimization variables are listed below.

Takes in PGM and applies mask to Flags
   - geometryMask
     Default value: 0.0

Limits number of geometry updates:
   - ItermaxPGM
     Default value: 0.0

Breaks simulation loop when SOR iterations fall below listed threshold:
   - sorIterationsThreshold
     Default value: 0.0

The following flags indicate the type of material modification allowed:
   - checkPressure
     Default value: 0.0

   - checkVelocity
     Default value: 0.0

   - checkVortex
     Default value: 0.0

The following flag determines if downstream surfaces are identified and penalized for material addition:

   - checkUpstream
     Default value: 0.0

Determines the factor by which downstream material addition is penalized:

   - downstreamVelocityFactor
     Default value: 0.0 # Warning please set this value if checkUpstream is set to 1

Pressure and velocity differences calculated across corner cells. Percentage threshold used to determine material removal:

   - percentPressure
     Default value: 0.0 # Warning please set this value if checkPressure is set to 1

   - percentVelocity
     Default value: 0.0 # Warning please set this value if checkVelocity is set to 1

Threshold values for addition and removal of cells based on velocity:

   - minVelocity
     Default value: 0.0 # Set to 0 if undesired

   - maxVelocity
     Default value: 0.0 # Set to large number if undesired

Fractional of the entire domain area which can be occupied by obstacles (setting=x --> area allowed=domain_size/x)

   - obstacleBudgetFraction 
     Default value: 0.0 # Warning please set this to an integer. Leaving this as 0 will break simulation

Minimum size of filled vortexes, vortexes smaller than given size (in grid cells), will be ignored

   - vortexSizeThreshold
     Default value: 0.0 # Highly recommended to set this if isVortex is 1. Detection of small/slow vortexes 
                            typically leads to instability.

Average velocity of filled vortex, vortexes with averaged velocities smaller than this value will be ignore.

   - vortexStrengthThreshold
     Default value: 0.0 # Highly recommended to set this if isVortex is 1. Detection of small/slow vortexes 
                            typically leads to instability.

                            
## How to use ./create-geometry.sh
This script converts JPEG files to the "compact" PGM format.  
As it is basically a wrapper around the "convert" command of the imagemagick suite, this will need to be installed
("imagemagick" package on Ubuntu).

The JPEG will need to be a canvas of the same size (WxH) of the destination PGM geometry, with colors:  

    black = obstacle  
    white = fluid  

The outer boundaries must not be represented in the JPEG, so its size must be equal to that of the inner domain.

The output file will have the same path and name of the input one, but will have the .pgm extension. Be careful not to
overwrite existing files!

The basic usage is:  
    `$ ./create-geometry.sh Scenario/geometry.jpg`  
which will create the Scenario/geometry.pgm file.

An experimental rescaling feature is also available, but be aware that it might introduce forbidden geometry
configurations! The following command will rescale to a 200x100 domain:  
    `$ ./create-geometry.sh Scenario/geometry.jpg 200 100`

The tool can also be used to convert a PGM back to the JPEG format, just by feeding the PGM file:  
    `$ ./create-geometry.sh Scenario/geometry.pgm`

