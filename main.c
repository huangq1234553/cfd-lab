#include "helper.h"
#include "visual.h"
#include "init.h"
#include "sor.h"
#include "boundary_val.h"
#include "uvp.h"
// #include "logger.h"
#include <mpi.h>
#include "parallel.h"
#include "logger.h"


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

int main(int argc, char **argv)
{
    
    // Handling the problem file name which is passed as 1st argument.
    char szFileName[256]; // We assume name will not be longer than 256 chars...
    strcpy(szFileName, argv[1]);
    strcat(szFileName, ".dat");
    //
    char problem[256];
    char geometry[1024]; // bigger since this can be a full path
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
    double eps;               /* accuracy bound for pressure*/
    double dt_value;          /* time for output */
    int n = 0;                  /* timestep iteration counter */
    double res = 10;          /* residual */
    double res_global;            /* used to store global residual*/
    double t = 0;              /* initial time */
    int it;                      /* sor iteration counter */
    int iproc;                  /* number of processes in i direction*/
    int jproc;                  /* number of processes in j direction*/
    double mindt = 10000;
    
    read_parameters(szFileName, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax,
                    &alpha, &omg,
                    &tau, &itermax, &eps, &dt_value, problem, geometry, &iproc, &jproc);
    strcpy(problem, argv[1]); // Just because we are not currently reading it from config. TODO: fix this properly.
    
    MPI_Status status;
    int my_rank, num_proc;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    int il, ir, jb, jt;
    int rank_l, rank_r, rank_b, rank_t;
    int omg_i, omg_j;
    
    // Now check if the declared processor grid (iproc x jproc) and the MPI num processes agree
    if (iproc*jproc != num_proc)
    {
        // TODO: here log a message also in the log/stdout
        ERROR("Mismatch between iproc*jproc and MPI number of processes.\n"
              "Please review your configuration and command-line to make sure they agree!");
        return 1;
    }
    
    init_parallel(iproc, jproc, imax, jmax, my_rank, &il, &ir, &jb, &jt, &rank_l, &rank_r, &rank_b, &rank_t,
                  &omg_i, &omg_j, num_proc);

    int imax_local = ir - il, jmax_local = jt - jb;

    if (my_rank == 0) {
        setLoggerOutputFolder("Out");
        setLoggerStartTime();
        openLogFile(); // Initialize the log file descriptor.
    }
    
    double **P = matrix(0, imax_local + 2, 0, jmax_local + 2);
    double **U = matrix(0, imax_local + 3, 0, jmax_local + 2);
    double **V = matrix(0, imax_local + 2, 0, jmax_local + 3);
    double **F = matrix(0, imax_local + 3, 0, jmax_local + 2);
    double **G = matrix(0, imax_local + 2, 0, jmax_local + 3);
    double **RS = matrix(0, imax_local, 0, jmax_local);
    double *bufSend = (double *) malloc((size_t) (3 * (imax_local + 3) * sizeof(double *)));
    double *bufRecv = (double *) malloc((size_t) (3 * (imax_local + 3) * sizeof(double *)));
    
    init_uvp(UI, VI, PI, imax_local, jmax_local, U, V, P);

    // Todo: This shouldn't be necessary and can probably be removed now
    // Perform one dt calculation prior to parallelized dt calculation in while loop.
    if (tau > 0)
    {
        dt = fmin(fmin((Re / 2 / (1 / pow(dx, 2) + 1 / pow(dy, 2))), fmin(dx / UI, dy / VI)), dt_value);
    }

//    printf("[R%d] Right before writing the viz...\n", my_rank); //debug
    if(my_rank == 0){
        logEvent(INFO, t, "Writing visualization file n = %d", n);
        //printf("INFO: Writing visualization file n=%d\n", n);
    }
    write_vtkFile(problem, n, my_rank, xlength, ylength, il, jb, imax_local+1, jmax_local+1, dx, dy, U, V, P);
    n++;
//    printf("[R%d] Right after writing the viz...\n", my_rank); //debug
    // Max values for U and V, for dt calculation
    double uMax = UI;
    double vMax = VI;
// 	// simulation interval 0 to t_end
    double currentOutputTime = 0; // For chosing when to output
    while (t < t_end)
    {
        // ensure boundary conditions for velocity
        // here we only need to set boundary values if local boundaries coincide with global boundaries
        if (omg_i == 0 || omg_i == iproc - 1 || omg_j == 0 || omg_j == jproc - 1)
        {
            boundaryvalues(omg_i, omg_j, iproc, jproc, imax_local, jmax_local, U, V);
        }
        
        // momentum equations M1 and M2 - F and G are the terms arising from explicit Euler velocity update scheme
        calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax_local, jmax_local, U, V, F, G);
        
        if (omg_i == 0 || omg_i == iproc - 1 || omg_j == 0 || omg_j == jproc - 1)
        {
            boundaryvalues_FG(omg_i, omg_j, iproc, jproc, imax_local, jmax_local, F, G, U, V);
        }

// 		// momentum equations M1 and M2 are plugged into continuity equation C to produce PPE - depends on F and G - RS is the rhs of the implicit pressure update scheme
        calculate_rs(dt, dx, dy, imax_local, jmax_local, F, G, RS);

// 		// solve the system of eqs arising from implicit pressure uptate scheme using succesive overrelaxation solver
        it = 0;
        res = 1e9;
        res_global = 1e9;
        while (it < itermax && res_global > eps)
        {
    
            sor(omg, dx, dy, imax_local + 1, jmax_local + 1, P, RS, &res);
            pressure_comm(P, rank_l, rank_r, rank_b, rank_t, bufSend, bufRecv, &status, imax_local + 1, jmax_local + 1);
            if (omg_i == 0 || omg_i == iproc - 1 || omg_j == 0 || omg_j == jproc - 1)
            {
                boundaryvalues_P(omg_i, omg_j, iproc, jproc,  imax_local + 1, jmax_local + 1, P);
            }
            it++;
            MPI_Allreduce(&res, &res_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            res_global = res_global/(imax*jmax);
            res_global = sqrt(res_global);
        }
        if (my_rank == 0)
        {
            logEvent(INFO, t, (char*)"dt=%f, numSorIterations=%d, sorResidual=%f", dt, it, res_global);
            if (it == itermax)
            {
                logEvent(WARNING, t, "Max number of iterations reached on SOR. Probably it did not converge!");
            }
        }

// 		// calculate velocities acc to explicit Euler velocity update scheme - depends on F, G and P
        calculate_uv(dt, dx, dy, imax_local, jmax_local, omg_i, omg_j, iproc, jproc, U, V, F, G, P, &uMax,
                     &vMax); // Here we get the uMax and vMax
        uv_comm(U, V, rank_l, rank_r, rank_b, rank_t, bufSend, bufRecv, &status, imax_local, jmax_local);
        
        // write visualization file for current iteration (only every dt_value step)
        if (t >= currentOutputTime)
        {
//            logEvent(t, "INFO: Writing visualization file n=%d", n);
            if(my_rank == 0){
                logEvent(INFO, t, "Writing visualization file n = %d", n);
                //printf("INFO: Writing visualization file n=%d\n", n);
            }
            write_vtkFile(problem, n, my_rank, xlength, ylength, il, jb, imax_local+1, jmax_local+1, dx, dy, U, V, P);
            currentOutputTime += dt_value;
            // update output timestep iteration counter
            n++;
        }
        // Recap shell output
        //        printf("[%12.9f] INFO: dt=%f, numSorIterations=%d, sorResidual=%f\n", t, dt, it, res);
//        logEvent(t, "INFO: dt=%f, numSorIterations=%d, sorResidual=%f", dt, it, res);

        // adaptive stepsize control based on stability conditions ensures stability of the method!
        // dt = tau * min(cond1, cond2, cond3) where tau is a safety factor
        // NOTE: if tau<0, stepsize is not adaptively computed!
        if (tau > 0)
        {
            double dtLocal;
            calculate_dt(Re, tau, &dtLocal, dx, dy, uMax, vMax);
            dtLocal = fmin(dtLocal, dt_value); // test, to avoid a dt bigger than visualization interval
//            printf("[R%d] Local dt before allreduce: %f\n", my_rank, dtLocal); //debug
            MPI_Allreduce(&dtLocal, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
//            printf("[R%d] Global dt after allreduce: %f\n", my_rank, dt); //debug
            // Used to check the minimum time-step for convergence
            if (dt < mindt)
            {
                mindt = dt;
            }
        }
        // advance in time
        t += dt;
    }
    
    // write visualisation file for the last iteration
    if (my_rank == 0) {
        logEvent(INFO, t, (char*)"Writing visualization file n = %d", n);
    }
    write_vtkFile(problem, n, my_rank, xlength, ylength, il, jb, imax_local+1, jmax_local+1, dx, dy, U, V, P);

    if (my_rank == 0) {
        // Check value of U[imax/2][7*jmax/8] (task6)
        //printf("Final value for U[imax/2][7*jmax/8] = %16e\n", U[imax / 2][7 * jmax / 8]);
        //logMsg("Final value for U[imax/2][7*jmax/8] = %16e", U[imax / 2][7 * jmax / 8]);
        printf("Min dt value used: %16e\n", mindt);
        logMsg("Min dt value used: %16e", mindt);
        closeLogFile(); // Properly close the log file
    }
    
    // Check value of U[imax/2][7*jmax/8] (task6)
    // logMsg("Final value for U[imax/2][7*jmax/8] = %16e", U[imax / 2][7 * jmax / 8]);
    free_matrix(P, 0, imax_local + 2, 0, jmax_local + 2);
    free_matrix(U, 0, imax_local + 3, 0, jmax_local + 2);
    free_matrix(V, 0, imax_local + 2, 0, jmax_local + 3);
    free_matrix(F, 0, imax_local + 3, 0, jmax_local + 2);
    free_matrix(G, 0, imax_local + 2, 0, jmax_local + 3);	


    free_matrix(RS, 0, imax_local, 0, jmax_local);
    free(bufSend);
    free(bufRecv);

    MPI_Finalize();

    
    return 0;
}

