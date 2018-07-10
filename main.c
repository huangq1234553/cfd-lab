#include <libgen.h>
#include <sys/stat.h>
#include <stdbool.h>
#include "helper.h"
#include "visual.h"
#include "init.h"
#include "sor.h"
#include "uvp.h"
#include "logger.h"
#include "timing.h"
#include <time.h>
#include <stdlib.h>

/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_flag(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the 
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */

double performSimulation(const char *outputFolder, const char *outputFolderPGM, const char *problem, double Re, double GX,
                         double GY, double t_end, double xlength, double ylength, double dt, double dx, double dy, int imax,
                         int jmax, double alpha, double omg, double tau, int itermax, double eps, double dt_value, int n,
                         int k, double res, double t, int it, double mindt, int noFluidCells, double beta, double Pr,
                         BoundaryInfo *boundaryInfo, double dt_check, int **Flags, double **U, double **V, double **F,
                         double **G, double **RS, double **P, double **T, int **PGM, bool computeTemperatureSwitch,
                         int itThreshold, double *maxU, double *maxV, double t_threshold, int *globalVtkCounter);

int main(int argc, char **argv)
{
    const char* USAGE = "\nUSAGE:\t./sim configuration.dat [OPTIONS]\n"
                        "\nOPTIONS:\n"
                        "\t-o PATH\t\tSet the output folder path (default is INPUT_DIR/Out).\n"
                        "\t--q\t\tSet the logging level to PRODUCTION (ERROR, WARNING and PRODUCTION traces will be enabled).\n"
                        "\t--debug\t\tSet the logging level to DEBUG (all traces will be enabled).\n"
                        "\t--notemp\tDisable temperature computation.\n"
                        "\t--fix-initial-geometry\n\t\t\tAllows auto-fixing the initial geometry before starting the simulation.\n"
                        "\n"
                        "INPUT_DIR is the directory where the configuration.dat file is placed.\n"
                        "Default logging level is INFO.\n"
                        "Logging levels ordered by decreasing priority are ERROR, WARNING, PRODUCTION, INFO, DEBUG.";
    // As first thing, initialize logger start time
    setLoggerStartTime();
    
    // Handling the problem file name which is passed as 1st argument.
    char szFileName[256] = ""; // We assume name will not be longer than 256 chars...
    RunningMode runningMode = EXTENDED;
    char computeTemperatureSwitch = 1;
    char inputFolder[512] = "";
    char outputFolder[512] = ""; // Please don't use superlong paths here, 512 should be more than enough :P
    char outputFolderPGM[512] = "";
    bool fixInitialGeometryAllowed = false;
    bool pikachuMode = false; //Just for fun
    
    //
    // BEGIN command line flags management
    //
    int requiredArgCount = 0;
    int i = 1; // Arg counter
    while (i < argc)
    {
        if (strcmp(argv[i], "--compact") == 0)
        {
            runningMode = COMPACT;
        }
        else if (strcmp(argv[i], "--notemp") == 0)
        {
            computeTemperatureSwitch = 0;
        }
        else if (strcmp(argv[i], "-o") == 0)
        {
            // Setting the output folder -- Here reading the following arg (so format is -o /path/to/folder )
            strcpy(outputFolder, argv[i + 1]);
            // Remove trailing slash if present
            if (outputFolder[strlen(outputFolder) - 1] == '/')
            {
                outputFolder[strlen(outputFolder) - 1] = 0;
            }
            ++i; // Extra increment since we read 2 arguments here
        }
        else if (strcmp(argv[i], "-q") == 0)
        {
            // Set quiet mode: DebugLevel -> PRODUCTION
            setLoggerDebugLevel(PRODUCTION);
        }
        else if (strcmp(argv[i], "--debug") == 0)
        {
            // Set debug mode: DebugLevel -> DEBUG
            setLoggerDebugLevel(DEBUG);
        }
        else if (strcmp(argv[i], "--fix-initial-geometry") == 0)
        {
            fixInitialGeometryAllowed = true;
        }
        else if (strcmp(argv[i], "--pikachu-mode") == 0)
        {
            // Just for fun :)
            pikachuMode = true;
            // Initialize random seed for random removal
            srand((unsigned int) time(NULL));
        }
        // All the below is just error catching... (new options must be set above here!)
        else if (argv[i][0] == '-')
        {
            char buf[1024];
            sprintf(buf, "Unrecognized option: %s\n%s", argv[i], USAGE);
            THROW_ERROR(buf);
        }
        else if (strlen(szFileName) != 0)
        {
            char buf[1024];
            sprintf(buf, "Unrecognized argument: %s\n%s", argv[i], USAGE);
            THROW_ERROR(buf);
        }
        else
        {
            strcpy(szFileName, argv[i]);
            // Now if argument doesn't end with ".dat", append extension to it
            size_t fnamelen = strlen(szFileName);
            if (fnamelen > 4 && strcmp(szFileName + fnamelen - 4, ".dat") != 0)
            {
                strcat(szFileName, ".dat");
            }
            // Extract folder to be used as input folder
            strcpy(inputFolder, szFileName);
            dirname(inputFolder);
            ++requiredArgCount;
        }
        ++i;
    }
    //
    if (requiredArgCount == 0)
    {
        char buf[1024];
        sprintf(buf, "No arguments passed!\n%s", USAGE);
        THROW_ERROR(buf);
    }
    
    // In case no outputFolder is passed, default to inputFolder/Out...
    if (strlen(outputFolder) == 0)
    {
        sprintf(outputFolder, "%s/Out", inputFolder);
    }
    // ...and create it if not on filesystem, in any case
    {
        struct stat st = {0};
        if (stat(outputFolder, &st) == -1)
        {
            mkdir(outputFolder, 0700);
        }
    }

    //always create a subfolder for PGMs outputFolder/PGM
    if (strlen(outputFolderPGM) == 0)
    {
        sprintf(outputFolderPGM, "%s/PGM", outputFolder);
    }
    // ...and create it if not on filesystem
    {
        struct stat st = {0};
        if (stat(outputFolderPGM, &st) == -1)
        {
            mkdir(outputFolderPGM, 0700);
        }
    }

    //
    // END command line flags management
    //
    setLoggerOutputFolder(outputFolder);
    openLogFile(); // Initialize the log file descriptor.
    //
    char problem[256];
    char geometry[512]; // bigger since this can be a full path
    char geometryMaskFile[512]; // bigger since this can be a full path
    double Re;                /* reynolds number   */
    double UI;                /* velocity x-direction */
    double VI;                /* velocity y-direction */
    double PI;                /* pressure */
    double GX;                /* gravitation x-direction */
    double GY;                /* gravitation y-direction */
    double t_end;             /* end time */
    double xlength;           /* length of the domain x-dir.*/
    double ylength;           /* length of the domain y-dir.*/
    double dt;                /* time step */
    double dx;                /* length of a cell x-dir. */
    double dy;                /* length of a cell y-dir. */
    int imax;                /* number of cells x-direction*/
    int jmax;                /* number of cells y-direction*/
    double alpha;             /* uppwind differencing factor*/
    double omg;               /* relaxation factor */
    double tau;               /* safety factor for time step*/
    int itermax;              /* max. number of iterations  */
    int itermaxPGM=50;           /* max. nunmber of iterations for PGM convergence */
    double eps;               /* accuracy bound for pressure*/
    double dt_value;          /* time for output */
    int n = 0;                  /* timestep iteration counter */
    int k = 1;                     /* PGM iteration counter */
    double res = 10;          /* residual */
    double t = 0;              /* initial time */
    int it = 0;                      /* sor iteration counter */
    double mindt = 10000;       /* arbitrary counter that keeps track of minimum dt value in calculation */
    int noFluidCells;          /* number of fluid cells in simulation */
    int noCouplingCells;       /* number of coupling cells in simulation */
    double beta;              /* coefficient of thermal expansion */
    double TI;                  /* initial temperature */
    double T_h;                  /* hot surface boundary condition */
    double T_c;                  /* cold surface boundary condition */
    double Pr;                  /* Prandtl number */
    double outflow = 0;
    double outflowMax = 0;
    double percentPressure = 0.6;
    double percentVelocity = 0.6;
    double minVelocity = 0.05;
    int isVortex, isPressure, isVelocity;
    double maxVelocity = 0.5;
    int obstacleBudget = 0;
    double obstacleBudgetFraction = 2;
    int isGeometryAdaptivityEnabled;
    int vortexAreaThreshold;
    double vortexStrengthThreshold;
    int isUpstreamCheckEnabled = 0;
    double downstreamVelocityFactor = 0;
    
    int maxK = 0;
    int sorIterationsAcceptanceThreshold;
    // Params for preCICE coupling
    double x_origin = 0.0, y_origin = 0.0;
    char precice_config[512], participant_name[128], mesh_name[128],
            read_data_name[128], write_data_name[128];
    
    BoundaryInfo boundaryInfo[4];
    
    read_parameters(szFileName, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax,
                    &alpha, &omg,
                    &tau, &itermax, &itermaxPGM, &sorIterationsAcceptanceThreshold, &eps, &dt_value, problem, geometry,
                    geometryMaskFile, boundaryInfo, &beta, &TI, &T_h, &T_c, &Pr, &x_origin, &y_origin,
                    &minVelocity, &downstreamVelocityFactor, &maxVelocity, &percentPressure, &percentVelocity, &isPressure, &isVelocity,
                    &isUpstreamCheckEnabled, &isVortex, &vortexAreaThreshold, &vortexStrengthThreshold, &obstacleBudgetFraction,
                    precice_config, participant_name, mesh_name, read_data_name, write_data_name, &isGeometryAdaptivityEnabled);
    
    // In case geometry was given as a filename only, prepend it with inputFolder path, just in case it is not CWD.
    if (strstr(geometry, "/") == NULL)
    {
        char buf[512];
        sprintf((char *) buf, "%s/%s", inputFolder, geometry);
        strcpy(geometry, buf);
        logMsg(PRODUCTION, "Using geometry file: %s", geometry);
    }
    // In case geometryMaskFile was given as a filename only, prepend it with inputFolder path, just in case it is not CWD.
    logMsg(DEBUG, "Raw GeometryMask filename: \"%s\"", geometryMaskFile);
    if (strstr(geometryMaskFile, "/") == NULL && !isStringDefault(geometryMaskFile))
    {
        char buf[512];
        sprintf((char *) buf, "%s/%s", inputFolder, geometryMaskFile);
        strcpy(geometryMaskFile, buf);
        logMsg(PRODUCTION, "Using geometryMask file: %s", geometryMaskFile);
    }
    
    double dt_check = fmin(dt, dt_value);
    
    int **Flags = imatrix(0, imax + 1, 0, jmax + 1);
    double **U = matrix(0, imax + 1, 0, jmax + 1);
    double **V = matrix(0, imax + 1, 0, jmax + 1);
    double **F = matrix(0, imax + 1, 0, jmax + 1);
    double **G = matrix(0, imax + 1, 0, jmax + 1);
    double **RS = matrix(0, imax + 1, 0, jmax + 1);
    double **P = matrix(0, imax + 1, 0, jmax + 1);
    double **T = matrix(0, imax + 1, 0, jmax + 1);
    int **PGM = imatrix(0, imax + 1, 0, jmax + 1);
    
    // create flag array to determine boundary conditions
    if (runningMode == COMPACT)
    {
        THROW_ERROR("NO COMPACT MODE ANYMORE, SORRY! :(");
    }
    else
    {
        logMsg(PRODUCTION, "Running in extended mode");
        read_boundary_parameters_extended_mode(szFileName, boundaryInfo, dx, dy, imax, jmax, geometry);
    }
    
    // Log if temperature is being computed or not
    if (!computeTemperatureSwitch)
    {
        logMsg(PRODUCTION, "NoTemp mode: Temperature is not computed");
    }
    
    init_flag(problem, geometry, geometryMaskFile, imax, jmax, Flags, &noFluidCells, &noCouplingCells, runningMode, fixInitialGeometryAllowed);
    int noGeometryCells = imax*jmax - noFluidCells;
    obstacleBudget = min((int) round((imax * jmax) / obstacleBudgetFraction), imax * jmax) - noGeometryCells;

    // initialise velocities and pressure
    init_uvpt(UI, VI, PI, TI, imax, jmax, U, V, P, T, Flags);
    double maxU, maxV;
    
    // simulation interval 0 to t_end
    long simulationStartTime = getCurrentTimeMillis();
    bool vortexDetectionEnabled = (isVortex==1);
    bool vortexDetectionEnabledTmp = vortexDetectionEnabled;
    int globalVtkCounter = 0;
    if (isGeometryAdaptivityEnabled == 1)
    {
//        int prevNoFluidCells = (int) round(noFluidCells/(1-0.04)); // Let's have some forced time at the beginnging (2s if t_end=50)
        int prevNoFluidCells = (int) round(noFluidCells/(1-0.3)); // Let's have some forced time at the beginnging (15s if t_end=50)
        for (k; k <= itermaxPGM; k++)
        {
            // saving the *.pgm
            logEvent(PRODUCTION, t,
                     "Writing PGM file k=%d, executionTime=%.3fs, outflowValue=%.3f, maxOutflow=%.3f at iteration %d",
                     k,
                     getTimeSpentSeconds(simulationStartTime, getCurrentTimeMillis()),
                     outflow,
                     outflowMax,
                     maxK
            );
            decode_flags(imax, jmax, Flags, PGM);
            write_pgm(imax + 2, jmax + 2, PGM, outputFolderPGM, problem, k);
            
            double geometryChangeFactor = fabs(noFluidCells - prevNoFluidCells) / prevNoFluidCells;
            geometryChangeFactor = fmin(geometryChangeFactor, 1);
            geometryChangeFactor = fmax(geometryChangeFactor, 0.01); // minimum requirement to avoid too short runs
            logMsg(PRODUCTION, "Geometry change factor [k=%d]: %f", k, geometryChangeFactor);
            prevNoFluidCells = noFluidCells;
            mindt = performSimulation(outputFolder, outputFolderPGM, problem, Re, GX, GY, t_end, xlength, ylength, dt,
                                      dx, dy, imax, jmax, alpha, omg, tau, itermax, eps, dt_value, n, k, res, t, it,
                                      mindt, noFluidCells, beta, Pr, boundaryInfo,
                                      dt_check, Flags, U, V, F, G, RS, P, T, PGM, computeTemperatureSwitch,
                                      sorIterationsAcceptanceThreshold, &maxU, &maxV,
                                      round(t_end * geometryChangeFactor), &globalVtkCounter);
        
            //update PGM here - go through all the flags and decide what needs to be changed and what not
            if (vortexDetectionEnabled && vortexDetectionEnabledTmp && (obstacleBudget <= vortexAreaThreshold)) //debug additional condition
            {
                // Once we reach budget 0, disable vortex detection permanently, as it can introduce instability
                vortexDetectionEnabledTmp = 0;
                logMsg(PRODUCTION, "Vortex detection disabled. k=%d", k);
            }
            else if (vortexDetectionEnabled && !vortexDetectionEnabledTmp && (obstacleBudget > vortexAreaThreshold)) //debug additional condition
            {
                // Once we reach budget 0, disable vortex detection permanently, as it can introduce instability
                vortexDetectionEnabledTmp = 1;
                logMsg(PRODUCTION, "Vortex detection enabled. k=%d", k);
            }
            if (vortexDetectionEnabled && vortexDetectionEnabledTmp && k%10==0)
            {
                expandVortexSeeds(imax, jmax, &noFluidCells, U, V, P, Flags, vortexAreaThreshold,
                                  vortexStrengthThreshold, &obstacleBudget);
                geometryFix(U, V, P, Flags, imax, jmax, &noFluidCells, &obstacleBudget);
            }
            if (pikachuMode)
            {
                double removalProbability = 0; //default
                if (k > round(itermaxPGM*0.1)
                    && k <= round(itermaxPGM*0.25))
                {
                    removalProbability = 0.05; //do nothing at the beginning
                }
                else if (k >= round(itermaxPGM*0.25)
                         && k <= round(itermaxPGM*0.5))
                {
                    removalProbability = 0.1; //increase
                    minVelocity /= 2;
                }
                else if (k >= round(itermaxPGM*0.5)
                        && k <= round(itermaxPGM*0.7))
                {
                    removalProbability = 0.5; //increase
                    minVelocity /= 4;
                }
                else if (k >= round(itermaxPGM*0.7)
                         && k <= round(itermaxPGM*0.85))
                {
                    removalProbability = 0.80; //increase
                    minVelocity = 0;
                }
                else if (k >= round(itermaxPGM*0.85)
                         && k <= round(itermaxPGM*0.99))
                {
                    removalProbability = 0.95; //clear everything at the end
                    minVelocity = 0;
                }
                else if (k >= round(itermaxPGM*0.99))
                {
                    removalProbability = 1; //clear everything at the end
                    minVelocity = 0;
                }                logMsg(WARNING,"Pikachu mode enabled, removing some random solid cells! :) - removalProbability=%f, minVelocity=%f", removalProbability, minVelocity);
                randomGeometryRemoval(imax, jmax, &noFluidCells, &obstacleBudget, Flags, P, U, V, removalProbability);
            }
            
            update_pgm(imax, jmax, &noFluidCells, Flags, P, U, V, minVelocity, maxVelocity, percentPressure,
                       percentVelocity, maxU, maxV, k, &obstacleBudget, dx, dy, isPressure, isVelocity,
                       isUpstreamCheckEnabled, downstreamVelocityFactor);
            //fix forbidden geometry in case it exists
            geometryFix(U, V, P, Flags, imax, jmax, &noFluidCells, &obstacleBudget);
            //fix velocities on internal geometries, for nice vtks
            velocityFix(U, V, Flags, imax, jmax);
        
            outputCalculation(U, V, Flags, imax, jmax, &outflow);
        
            if (outflow > outflowMax)
            {
                outflowMax = outflow;
                maxK = k;
            }
        }
    }
    else
    {
        // In case of no geometry adaptivity, just run the simulation
        mindt = performSimulation(outputFolder, outputFolderPGM, problem, Re, GX, GY, t_end, xlength, ylength, dt, dx,
                                  dy, imax, jmax, alpha, omg, tau, itermax, eps, dt_value, n, k, res, t, it, mindt,
                                  noFluidCells, beta, Pr, boundaryInfo,
                                  dt_check, Flags, U, V, F, G, RS, P, T, PGM, computeTemperatureSwitch,
                                  sorIterationsAcceptanceThreshold, &maxU, &maxV, 0, &globalVtkCounter);
    }

    long simulationEndTime = getCurrentTimeMillis();

    // Check value of U[imax/2][7*jmax/8] (task6)
    logMsg(PRODUCTION, "Final value for U[imax/2][7*jmax/8] = %16e", U[imax / 2][7 * jmax / 8]);
    logMsg(PRODUCTION, "Min dt value used: %16e", mindt);
    logMsg(PRODUCTION, "Total time spent in simulation: %.3f s", getTimeSpentSeconds(simulationStartTime, simulationEndTime));
    
    free_imatrix(Flags, 0, imax + 1, 0, jmax + 1);
    free_matrix(U, 0, imax + 1, 0, jmax + 1);
    free_matrix(V, 0, imax + 1, 0, jmax + 1);
    free_matrix(F, 0, imax + 1, 0, jmax + 1);
    free_matrix(G, 0, imax + 1, 0, jmax + 1);
    free_matrix(RS, 0, imax + 1, 0, jmax + 1);
    free_matrix(P, 0, imax + 1, 0, jmax + 1);
    free_matrix(T, 0, imax + 1, 0, jmax + 1);
    freeAllBoundaryInfo(boundaryInfo);
    
    closeLogFile(); // Properly close the log file
    
    return 0;
}

double performSimulation(const char *outputFolder, const char *outputFolderPGM, const char *problem, double Re, double GX,
                         double GY, double t_end, double xlength, double ylength, double dt, double dx, double dy, int imax,
                         int jmax, double alpha, double omg, double tau, int itermax, double eps, double dt_value, int n,
                         int k, double res, double t, int it, double mindt, int noFluidCells, double beta, double Pr,
                         BoundaryInfo *boundaryInfo, double dt_check, int **Flags, double **U, double **V, double **F,
                         double **G, double **RS, double **P, double **T, int **PGM, bool computeTemperatureSwitch,
                         int itThreshold, double *maxU, double *maxV, double t_threshold, int *globalVtkCounter)
{
    double currentOutputTime = 0; // For chosing when to output
    long interVisualizationExecTimeStart = getCurrentTimeMillis();
    it = itThreshold + 1;
    setPressureOuterBoundaryValues(imax, jmax, P, Flags, boundaryInfo);
    while ((it > itThreshold || t < t_threshold) && t < t_end)
    {
        
        // adaptive stepsize control based on stability conditions ensures stability of the method!
        // dt = tau * min(cond1, cond2, cond3) where tau is a safety factor
        // NOTE: if tau<0, stepsize is not adaptively computed!
        if (tau > 0)
        {
            calculate_dt(Re, Pr, tau, &dt, dx, dy, imax, jmax, U, V, maxU, maxV);
            dt = fmin(dt, dt_check); // test, to avoid a dt bigger than visualization interval
            // Used to check the minimum time-step for convergence
            if (dt < mindt)
            {
                mindt = dt;
            }
        }
        
        // ensure boundary conditions for velocity
        // Special boundary condition are addressed here by using the boundaryInfo data.
        // These special boundary values are configured at configuration time in read_parameters().
        boundaryval(imax, jmax, U, V, T, Flags, boundaryInfo);
        
        // calculate T using energy equation in 2D with boussinesq approximation
        if (beta != 0 &&
            computeTemperatureSwitch)
        { // If beta==0 we won't add the temperature term in F and G, so we don't need to compute it!
            calculate_T(Re, Pr, dt, dx, dy, alpha, imax, jmax, T, U, V, Flags);
        }
        
        // momentum equations M1 and M2 - F and G are the terms arising from explicit Euler velocity update scheme
        calculate_fg(Re, GX, GY, alpha, beta, dt, dx, dy, imax, jmax, U, V, F, G, T, Flags);
        
        // momentum equations M1 and M2 are plugged into continuity equation C to produce PPE - depends on F and G - RS is the rhs of the implicit pressure update scheme
        calculate_rs(dt, dx, dy, imax, jmax, F, G, RS, Flags);
        
        // solve the system of eqs arising from implicit pressure uptate scheme using succesive overrelaxation solver
        it = 0;
        res = 1e9;
        while (it < itermax && res > eps)
        {
            sor(omg, dx, dy, imax, jmax, P, RS, Flags, boundaryInfo, &res, noFluidCells);
            it++;
        }
        // calculate velocities acc to explicit Euler velocity update scheme - depends on F, G and P
        calculate_uv(dt, dx, dy, imax, jmax, U, V, F, G, P, Flags);
        
        // write visualization file for current iteration (only every dt_value step)
        if (t >= currentOutputTime)
        {
            logEvent(PRODUCTION, t, "Writing visualization file k=%d, n=%d, executionTime=%.3fs",
                     k, n,
                     getTimeSpentSeconds(interVisualizationExecTimeStart, getCurrentTimeMillis())
            );
            write_vtkFile(outputFolder, problem, n, k, xlength, ylength, imax, jmax, dx, dy, U, V, P, T, Flags);
            if (k==1)
            {
                // Now also add a vtk to a series of global streams (PGM id = -1)
                write_vtkFile(outputFolder, problem, (*globalVtkCounter)++, -1, xlength, ylength, imax, jmax, dx, dy, U,
                              V, P, T, Flags);
            }
    
    
            currentOutputTime += dt_value;
            // update output timestep iteration counter
            n++;
            interVisualizationExecTimeStart = getCurrentTimeMillis();
        }
        // Recap shell output
        if (it == itermax)
        {
            logEvent(WARNING, t, "Max number of iterations reached on SOR. Probably it did not converge!");
            logEvent(WARNING, t, "dt=%f, numSorIterations=%d, sorResidual=%f, noFluidCells=%d", dt, it, res, noFluidCells);
            
        }
        else
        {
            logEvent(INFO, t, "dt=%f, numSorIterations=%d, sorResidual=%f, noFluidCells=%d", dt, it, res, noFluidCells);
        }
        // advance in time
        t += dt;
    }
    
    // write visualisation file for the last iteration
    logEvent(PRODUCTION, t, "Writing visualization file k=%d, n=%d, executionTime=%.3fs",
             k, n,
             getTimeSpentSeconds(interVisualizationExecTimeStart, getCurrentTimeMillis())
    );
    write_vtkFile(outputFolder, problem, n, k, xlength, ylength, imax, jmax, dx, dy, U, V, P, T, Flags);
    // Now also add a vtk to a series of "asymptotic" streams (PGM id = 0)
    write_vtkFile(outputFolder, problem, k, 0, xlength, ylength, imax, jmax, dx, dy, U, V, P, T, Flags);
    // Now also add a vtk to a series of global streams (PGM id = -1)
    write_vtkFile(outputFolder, problem, (*globalVtkCounter)++, -1, xlength, ylength, imax, jmax, dx, dy, U, V, P,
                  T, Flags);

    return mindt;
}

//eof
