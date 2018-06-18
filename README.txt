### README - How to execute the code and the test cases ###

1) Build:

    $ make

2) Execute the predefined simulation cases via make commands:
	The make commands for the various cases are:
	solidplate: 		$ make runsolidplate
	solidconvection: 	$ make runsolidconvection
	solidexchange: 		$ make runsolidexchangeF1 	& 	$ make runsolidexchangeF2 
	
    Remember to execute the respective OpenFOAM solver in another terminal:
    solidplate:         $ ./run_solid_plate
    solidconvection:    $ ./run_solid_convection
    solidexchange:      $ ./run_solid_exchange

The .vtk files are output into the Solid_*/Out* folders.


### Supported command line arguments (OPTIONS)
List of supported arguments:
    --notemp
        This disables the temperature computations, in case they are not required.

    -q
        "Quiet" mode: reduces the verbosity of output/log by disabling INFO level traces (it still retains PRODUCTION,
        WARNING, ERROR level traces).

    -o path/to/out/folder/
        Output directory path. In case this is not set, an "Out" subfolder will be created inside the folder containing
        the configuration .dat file.

### Logging
The same output which is written to stdout is also written into a log file (sim.log) in the output folder, along with
the output visualization files.
This is useful to keep trace of what parameters where used to generate the output and to investigate potential issues.

### Configuring boundary parameters (for velocity and temperature) in the .dat configuration file
Setting velocity boundary parameters and temperature boundary parameters
requires setting extra variables in configuration file.
If a variable is omitted in the configuration file, it will take its default value.

    - left_boundary_dirichlet_U, right_boundary_dirichlet_U, top_boundary_dirichlet_U, bottom_boundary_dirichlet_U
    - left_boundary_dirichlet_V, right_boundary_dirichlet_V, top_boundary_dirichlet_V, bottom_boundary_dirichlet_V
            Accepted values: any double
            Default value: 0.0
            
    - left_boundary_T, right_boundary_T, top_boundary_T, bottom_boundary_T
    - left_boundary_qN, right_boundary_qN, top_boundary_qN, bottom_boundary_qN
            Accepted values: any double
            Default value: 0.0

    - left_boundary_k, right_boundary_k, top_boundary_k, bottom_boundary_k
                Accepted values: any double
                Default value: 1.0
#eof
