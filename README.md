# cfd-lab

Repository for the code of GroupA, CFD Lab course SoSe2018 at TUM.
Enjoy! :)

## README - How to execute the code and the test cases

1) Easy way to get started: execute the pre-defined test cases!

    `$ make`  
    `$ make tests`

2) Configuring and executing your own scenario:

    `$ mkdir Scenario`  
    `$ touch Scenario/scenario.dat`

    Now edit the scenario.dat file to contain all the configuration for the scenario (more details below).  
    Then create a geometry .pgm file (below instruction for the ./create-geometry.sh script) and refer to it in the scenario.dat.

    (Don't forget to build:
        `$ make`
    )

    Then the scenario can be executed with:  
    `$ ./sim Scenario/scenario.dat [OPTIONS]`

## Supported command line arguments (OPTIONS)
List of supported arguments:

    --compact
        This is required in order to use "compact" PGM geometry files (more info on this below).
        If not set, the program will assume an "extended" PGM file.

    --notemp
        This disables the temperature computations, in case they are not required.

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

## Supported PGM "formats"
Our code supports 2 types of PGM formats:
1) Our own "compact" one, where:
  - the PGM only represents the inner domain (so outer boundary halo is left out)
  - in the PGM a value of 1 is an obstacle cell and 0 is a fluid cell

2) The "extended" format as per example published on Moodle, where:
  - the domain outer halo is included in the PGM file
  - 0=NOSLIP, 1=FREESLIP, 2=OUTFLOW, 3=INFLOW, 4=FLUID

## Configuring boundary parameters (for velocity and temperature) in the .dat configuration file
Setting velocity boundary parameters (when running in "compact" mode) and temperature boundary parameters (always)
requires setting extra variables in configuration file.
If a variable is omitted in the configuration file, it will take its default value.

  - left_boundary_type, right_boundary_type, top_boundary_type, bottom_boundary_type
    Accepted values: NOSLIP, FREESLIP, MOVINGWALL, INFLOW, OUTFLOW  
    Default value: NOSLIP

  - left_boundary_U, right_boundary_U, top_boundary_U, bottom_boundary_U  
    left_boundary_V, right_boundary_V, top_boundary_V, bottom_boundary_V  
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

