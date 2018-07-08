#ifndef __HELPER_H__
#define __HELPER_H__

/* includefiles */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <time.h>


#ifdef PI
#undef PI
#endif

#define FREE_ARG char*

/**
 * Maximum length of input lines 
 */
#define MAX_LINE_LENGTH 1024	   

/**
 * Stores the last timer value 
 */
//typedef enum Direction {CENTER=1, TOP=16, BOT=8, LEFT=4, RIGHT=2} Direction;
typedef enum Direction {CENTER=0, TOP=4, BOT=3, LEFT=2, RIGHT=1} Direction;
typedef enum FlagTypeBit {NSBIT=5, FSBIT=6, OFBIT=7, IFBIT=8, CBIT=9, TBIT=10, GEOMETRYMASKBIT=11} FlagTypeBit;
typedef enum FlippedReasonBit {GEOMETRY_FIX=12, VORTEX_FILL=13, } FlippedReasonBit;
typedef enum Optional {REQUIRED, OPTIONAL} Optional;
extern clock_t last_timer_reset;   


int min( int a, int b);	       
int max( int a, int b);
double fmin( double a, double b);
double fmax( double a, double b);
static inline short sign(int x)
{
    return (x > 0) - (x < 0);
}
static inline short fsign(double x)
{
    return (x > 0) - (x < 0);
}

int isObstacle(int flag);   // Current cell is an obstacle
int isFluid(int flag);      // Current cell is fluid
int isCoupling(int flag);      // Current cell is coupling
int isOutflow(int flag);      // Current cell is outflow
int isInflow(int flag);      // Current cell is inflow
int isGeometryConstant(int flag);      // Current cell does not allow to be switched from fluid to solid or viceversa
int isNeighbourObstacle(int flag, Direction direction); // Current cell's neighbor in the specified direction is obstacle
int isNeighbourFluid(int flag, Direction direction);    // Current cell's neighbor in the specified direction is fluid
int doesNeighbourAllowFlipToSolid(int **Flags, int i, int j, Direction direction, int imax, int jmax); // Current cell's neighbor in the specified direction can be flipped to solid
int hasAtLeastOneObstacleNeighbour(int flag);
int isCorner(int flag); // Current cell is a corner obstacle
int skipU(int flag);    // Current cell is surrounded by obstacles to its Top-Right-Bottom
int skipV(int flag);    // Current cell is surrounded by obstacles to its Left-Top-Right
void geometryCheck(int** flag, int imax, int jmax);  //Checks if forbidden geometry is in pgm
void decode_flags(int imax, int jmax, int **Flag, int** pic); // decode flags into pgm file gray scale

void write_pgm(int xsize, int ysize, int **pgm, const char *outputFolder, const char *szProblem, int iterationNumber); //write *.pgm file

void update_pgm(int imax, int jmax, int *noFluidCells, int **Flag, double **P, double **U, double **V, double minVelocity,
                double maxVelocity, double percentPressure, double percentVelocity, double maxU, double maxV, int k,
                int *obstacleBudget, double dx, double dy, int isPressure, int isVelocity, int isUpstreamCheckEnabled,
                double downstreamVelocityFactor);

void flipToFluid(double **U, double **V, int **Flags, int i, int j, int *obstacleBudget);
void flipToSolid(double **U, double **V, double **P, int **Flag, int i, int j, int *obstacleBudget);
double getRelativeVelocityThreshold(double maxU, double maxV, double percent);
int isVelocityGreaterThanRelativeThreshold(double velocity, double maxU, double maxV, double percent);
int isVelocityAboveRelativeThreshold(int isFlip, double percent, double **U, double **V, int Flag, int i, int j,
                                     double maxU, double maxV);
int isVelocityMagnitudeBelowThreshold(double eps, double leftVelocity, double rightVelocity, double bottomVelocity,
                                      double topVelocity);

int checkPressure(double percent, double ** P, int Flag, int i, int j);
void
geometryFix(double **U, double **V, double **P, int **Flag, int imax, int jmax, int *noFluidCells, int *obstacleBudget);
void velocityFix(double **U, double **V, int **Flag, int imax, int jmax);

void outputCalculation(double **U, double **V, int **Flags, int imax, int jmax, double *outflow);

int isGradientAtObstacleAboveThreshold(double **U, double **V, double dx, double dy, int cell, int i, int j,
                                       double gradThreshold);

void expandVortexSeeds(int imax, int jmax, int *noFluidCells, double **U, double **V, double **P, int **Flag,
                       int vortexAreaThreshold, double swirlThresholdStrength, int *obstacleBudget);

/**
 * Error handling:
 *
 * ERROR(s) writes an error message and terminates the program
 *
 * Example:
 * ERROR("File not found !");
 */
#define THROW_ERROR(s)    errhandler( __LINE__, __FILE__, s)

/**
 * Error handling:
 *
 * ERROR(s) writes an error message and terminates the program
 *
 * Example:
 * ERROR("File not found !");
 */
#define ERROUT stdout

/**
 * Error handling:
 *
 * ERROR(s) writes an error message and terminates the program
 *
 * Example:
 * ERROR("File not found !");
 */
void  errhandler( int nLine, const char *szFile, const char *szString );


/**
 * Reading from a datafile.
 *
 * The foloowing three macros help reading values from the parameter file.
 * If a variable cannot be found, the program stops with an error message.
 *
 * Example:
 * READ_INT( "MyFile.dat", imax );
 * READ_STRING( szFile, szProblem );
 */
#define READ_INT( szFileName, VarName, Optional)    read_int   ( szFileName, #VarName, &(VarName), Optional )

/**
 * Reading from a datafile.
 *
 * The foloowing three macros help reading values from the parameter file.
 * If a variable cannot be found, the program stops with an error message.
 *
 * Example:
 * READ_INT( "MyFile.dat", imax );
 * READ_STRING( szFile, szProblem );
 */
#define READ_DOUBLE( szFileName, VarName, Optional) read_double( szFileName, #VarName, &(VarName), Optional )

/**
 * Reading from a datafile.
 *
 * The foloowing three macros help reading values from the parameter file.
 * If a variable cannot be found, the program stops with an error message.
 *
 * Example:
 * READ_INT( "MyFile.dat", imax );
 * READ_STRING( szFile, szProblem );
 */
#define READ_STRING( szFileName, VarName, Optional) read_string( szFileName, #VarName,  (VarName), Optional )

void read_string(const char *szFilename, const char *szName, char *sValue, Optional optional);
void read_int(const char *szFilename, const char *szName, int *nValue, Optional optional);
void read_double(const char *szFilename, const char *szName, double *Value, Optional optional);


/**
 * Writing matrices to a file.
 * -----------------------------------------------------------------------
 * write_matrix(...) wites a matrice to a file
 * the file has the following format
 *
 *    -----------------------------------------
 *    |  xlength          |  float  |  ASCII  |
 *    ----------------------------------------|
 *    |  ylength          |  float  |  ASCII  |
 *    ----------------------------------------|
 *    |  nrl              |  int    |  ASCII  |
 *    ----------------------------------------|
 *    |  nrh              |  int    |  ASCII  |    1. call of the
 *    ----------------------------------------|
 *    |  ncl              |  int    |  ASCII  |
 *    ----------------------------------------|
 *    |  nch              |  int    |  ASCII  |    1. call of the
 *    ----------------------------------------|    function with
 *    |  m[nrl][ncl]      |  float  |  binaer |    bFirst == 1
 *    ----------------------------------------|
 *    |  m[nrl][ncl+1]    |  float  |  binaer |
 *    ----------------------------------------|
 *    |                  .                    |
 *                       .
 *    |                  .                    |
 *    -----------------------------------------
 *    |  m[nrh][nch]      |  float  |  binary |
 *    -----------------------------------------------------------------
 *    |  m[nrl][ncl]      |  float  |  binary |
 *    ----------------------------------------|
 *    |  m[nrl][ncl+1]    |  float  |  binary |     2. call with
 *    ----------------------------------------|     bFirst == 0
 *    |                  .                    |
 *                       .
 *    |                  .                    |
 *    -----------------------------------------
 *    |  m[nrh][nch]      |  float  |  binary |
 *    ------------------------------------------------------------------
 *
 * @param szFileName          name of the file
 * @param m                   matrix
 * @param nrl                 first column
 * @param nrh  		          last column
 * @param ncl                 first row
 * @param nch                 last row
 * @param xlength             size of the geometry in x-direction
 * @param ylength             size of the geometry in y-direction
 * @param xlength             size of the geometry in x-direction
 * @param fFirst              0 == append, else overwrite
 */
void write_matrix( 
  const char* szFileName,
  double **m,
  int nrl,
  int nrh,
  int ncl,
  int nch,
  double xlength,
  double ylength,	       
  int fFirst 
);

/**
 * @param szFileName    filehandle
 * @param m             matrix
 * @param nrl           first column
 * @param nrh           last column
 * @param ncl           first row
 * @param nch           last row
 */
void read_matrix( const char* szFileName,	               /* filehandle */
		  double **m,		       /* matrix */
		  int nrl,		       /* first column */
		  int nrh,		       /* last column */
		  int ncl,		       /* first row */
		  int nch );                   /* last row */


/**
 * matrix(...)        storage allocation for a matrix (nrl..nrh, ncl..nch)
 * free_matrix(...)   storage deallocation
 * init_matrix(...)   initialization of all matrix entries with a fixed
 *                  (floating point) value
 * imatrix(...)       analog for matrices with integer-entries
 *
 * Example:
 *    U = matrix ( 0 , imax+1 , 0 , jmax+1 );
 *    init_matrix( U , 0, imax+1, 0, jmax+1, 0 );
 *    free_matrix( U,  0, imax+1, 0, jmax+1 );
 */
double **matrix( int nrl, int nrh, int ncl, int nch );
/**
 * matrix(...)        storage allocation for a matrix (nrl..nrh, ncl..nch)
 * free_matrix(...)   storage deallocation
 * init_matrix(...)   initialization of all matrix entries with a fixed
 *                  (floating point) value
 * imatrix(...)       analog for matrices with integer-entries
 *
 * Example:
 *    U = matrix ( 0 , imax+1 , 0 , jmax+1 );
 *    init_matrix( U , 0, imax+1, 0, jmax+1, 0 );
 *    free_matrix( U,  0, imax+1, 0, jmax+1 );
 */
void free_matrix( double **m, int nrl, int nrh, int ncl, int nch );
/**
 * matrix(...)        storage allocation for a matrix (nrl..nrh, ncl..nch)
 * free_matrix(...)   storage deallocation
 * init_matrix(...)   initialization of all matrix entries with a fixed
 *                  (floating point) value
 * imatrix(...)       analog for matrices with integer-entries
 *
 * Example:
 *    U = matrix ( 0 , imax+1 , 0 , jmax+1 );
 *    init_matrix( U , 0, imax+1, 0, jmax+1, 0 );
 *    free_matrix( U,  0, imax+1, 0, jmax+1 );
 */
void init_matrix( double **m, int nrl, int nrh, int ncl, int nch, double a);

/**
 * matrix(...)        storage allocation for a matrix (nrl..nrh, ncl..nch)
 * free_matrix(...)   storage deallocation
 * init_matrix(...)   initialization of all matrix entries with a fixed
 *                  (floating point) value
 * imatrix(...)       analog for matrices with integer-entries
 *
 * Example:
 *    U = matrix ( 0 , imax+1 , 0 , jmax+1 );
 *    init_matrix( U , 0, imax+1, 0, jmax+1, 0 );
 *    free_matrix( U,  0, imax+1, 0, jmax+1 );
 */
int  **imatrix( int nrl, int nrh, int ncl, int nch );
/**
 * matrix(...)        storage allocation for a matrix (nrl..nrh, ncl..nch)
 * free_matrix(...)   storage deallocation
 * init_matrix(...)   initialization of all matrix entries with a fixed
 *                  (floating point) value
 * imatrix(...)       analog for matrices with integer-entries
 *
 * Example:
 *    U = matrix ( 0 , imax+1 , 0 , jmax+1 );
 *    init_matrix( U , 0, imax+1, 0, jmax+1, 0 );
 *    free_matrix( U,  0, imax+1, 0, jmax+1 );
 */
void free_imatrix( int **m, int nrl, int nrh, int ncl, int nch );
/**
 * matrix(...)        storage allocation for a matrix (nrl..nrh, ncl..nch)
 * free_matrix(...)   storage deallocation
 * init_matrix(...)   initialization of all matrix entries with a fixed
 *                  (floating point) value
 * imatrix(...)       analog for matrices with integer-entries
 *
 * Example:
 *    U = matrix ( 0 , imax+1 , 0 , jmax+1 );
 *    init_matrix( U , 0, imax+1, 0, jmax+1, 0 );
 *    free_matrix( U,  0, imax+1, 0, jmax+1 );
 */
void init_imatrix( int **m, int nrl, int nrh, int ncl, int nch, int a);


/**
 * reads in a ASCII pgm-file and returns the colour information in a two-dimensional integer array.
 */
int **read_pgm(const char *filename);


/**
 *                         useful macros
 * -----------------------------------------------------------------------
 *  The following macros can be helpful to display variables during the
 *  runtime of the program.
 *  If you start the program in a shell from xemacs, you can jump to the
 *  respectove rows by switching to the compilation-minor-mode.
 *
 *  DUMP_POSITION()           dumps the actual position within the program
 *  DUMP_MESSAGE( .)          dump a message in addition
 *  DUMP_INT(..)              dump an integer variable
 *
 *  DUMP_MATRIX_TO_FILE(..)
 *  DUMP_INT_TO_FILE(..)      writes the value of the variable in
 *  DUMP_DOUBLE_TO_FILE(..)   a tracefile
 *  DUMP_STRING_TO_FILE(..)
 *
 *  RESET_TIMER()     set timer to zero
 *  DUMP_TIMER()      dump time that has passed since the last
 *                    RESET_TIMER()
 */
#define DUMPOUT stdout

/**
 *                         useful macros
 * -----------------------------------------------------------------------
 *  The following macros can be helpful to display variables during the
 *  runtime of the program.
 *  If you start the program in a shell from xemacs, you can jump to the
 *  respectove rows by switching to the compilation-minor-mode.
 *
 *  DUMP_POSITION()           dumps the actual position within the program
 *  DUMP_MESSAGE( .)          dump a message in addition
 *  DUMP_INT(..)              dump an integer variable
 *
 *  DUMP_MATRIX_TO_FILE(..)
 *  DUMP_INT_TO_FILE(..)      writes the value of the variable in
 *  DUMP_DOUBLE_TO_FILE(..)   a tracefile
 *  DUMP_STRING_TO_FILE(..)
 *
 *  RESET_TIMER()     set timer to zero
 *  DUMP_TIMER()      dump time that has passed since the last
 *                    RESET_TIMER()
 */
#define DUMP_POSITION() fprintf( DUMPOUT, "%s:%d Dumpposition \n", __FILE__, __LINE__ )

/**
 *                         useful macros
 * -----------------------------------------------------------------------
 *  The following macros can be helpful to display variables during the
 *  runtime of the program.
 *  If you start the program in a shell from xemacs, you can jump to the
 *  respectove rows by switching to the compilation-minor-mode.
 *
 *  DUMP_POSITION()           dumps the actual position within the program
 *  DUMP_MESSAGE( .)          dump a message in addition
 *  DUMP_INT(..)              dump an integer variable
 *
 *  DUMP_MATRIX_TO_FILE(..)
 *  DUMP_INT_TO_FILE(..)      writes the value of the variable in
 *  DUMP_DOUBLE_TO_FILE(..)   a tracefile
 *  DUMP_STRING_TO_FILE(..)
 *
 *  RESET_TIMER()     set timer to zero
 *  DUMP_TIMER()      dump time that has passed since the last
 *                    RESET_TIMER()
 */
#define DUMP_MESSAGE(s) fprintf( DUMPOUT, "%s:%d %s\n",            __FILE__, __LINE__, s  )

/**
 *                         useful macros
 * -----------------------------------------------------------------------
 *  The following macros can be helpful to display variables during the
 *  runtime of the program.
 *  If you start the program in a shell from xemacs, you can jump to the
 *  respectove rows by switching to the compilation-minor-mode.
 *
 *  DUMP_POSITION()           dumps the actual position within the program
 *  DUMP_MESSAGE( .)          dump a message in addition
 *  DUMP_INT(..)              dump an integer variable
 *
 *  DUMP_MATRIX_TO_FILE(..)
 *  DUMP_INT_TO_FILE(..)      writes the value of the variable in
 *  DUMP_DOUBLE_TO_FILE(..)   a tracefile
 *  DUMP_STRING_TO_FILE(..)
 *
 *  RESET_TIMER()     set timer to zero
 *  DUMP_TIMER()      dump time that has passed since the last
 *                    RESET_TIMER()
 */
#define DUMP_INT(n)     fprintf( DUMPOUT, "%s:%d %s = %d\n", __FILE__, __LINE__, #n, n )

/**
 *                         useful macros
 * -----------------------------------------------------------------------
 *  The following macros can be helpful to display variables during the
 *  runtime of the program.
 *  If you start the program in a shell from xemacs, you can jump to the
 *  respectove rows by switching to the compilation-minor-mode.
 *
 *  DUMP_POSITION()           dumps the actual position within the program
 *  DUMP_MESSAGE( .)          dump a message in addition
 *  DUMP_INT(..)              dump an integer variable
 *
 *  DUMP_MATRIX_TO_FILE(..)
 *  DUMP_INT_TO_FILE(..)      writes the value of the variable in
 *  DUMP_DOUBLE_TO_FILE(..)   a tracefile
 *  DUMP_STRING_TO_FILE(..)
 *
 *  RESET_TIMER()     set timer to zero
 *  DUMP_TIMER()      dump time that has passed since the last
 *                    RESET_TIMER()
 */
#define DUMP_DOUBLE(d)  fprintf( DUMPOUT, "%s:%d %s = %f\n", __FILE__, __LINE__, #d, d )

/**
 *                         useful macros
 * -----------------------------------------------------------------------
 *  The following macros can be helpful to display variables during the
 *  runtime of the program.
 *  If you start the program in a shell from xemacs, you can jump to the
 *  respectove rows by switching to the compilation-minor-mode.
 *
 *  DUMP_POSITION()           dumps the actual position within the program
 *  DUMP_MESSAGE( .)          dump a message in addition
 *  DUMP_INT(..)              dump an integer variable
 *
 *  DUMP_MATRIX_TO_FILE(..)
 *  DUMP_INT_TO_FILE(..)      writes the value of the variable in
 *  DUMP_DOUBLE_TO_FILE(..)   a tracefile
 *  DUMP_STRING_TO_FILE(..)
 *
 *  RESET_TIMER()     set timer to zero
 *  DUMP_TIMER()      dump time that has passed since the last
 *                    RESET_TIMER()
 */
#define DUMP_STRING(s)  fprintf( DUMPOUT, "%s:%d %s = %s\n", __FILE__, __LINE__, #s, s )

/**
 *                         useful macros
 * -----------------------------------------------------------------------
 *  The following macros can be helpful to display variables during the
 *  runtime of the program.
 *  If you start the program in a shell from xemacs, you can jump to the
 *  respectove rows by switching to the compilation-minor-mode.
 *
 *  DUMP_POSITION()           dumps the actual position within the program
 *  DUMP_MESSAGE( .)          dump a message in addition
 *  DUMP_INT(..)              dump an integer variable
 *
 *  DUMP_MATRIX_TO_FILE(..)
 *  DUMP_INT_TO_FILE(..)      writes the value of the variable in
 *  DUMP_DOUBLE_TO_FILE(..)   a tracefile
 *  DUMP_STRING_TO_FILE(..)
 *
 *  RESET_TIMER()     set timer to zero
 *  DUMP_TIMER()      dump time that has passed since the last
 *                    RESET_TIMER()
 */
#define RESET_TIMER()   last_timer_reset = clock()

/**
 *                         useful macros
 * -----------------------------------------------------------------------
 *  The following macros can be helpful to display variables during the
 *  runtime of the program.
 *  If you start the program in a shell from xemacs, you can jump to the
 *  respectove rows by switching to the compilation-minor-mode.
 *
 *  DUMP_POSITION()           dumps the actual position within the program
 *  DUMP_MESSAGE( .)          dump a message in addition
 *  DUMP_INT(..)              dump an integer variable
 *
 *  DUMP_MATRIX_TO_FILE(..)
 *  DUMP_INT_TO_FILE(..)      writes the value of the variable in
 *  DUMP_DOUBLE_TO_FILE(..)   a tracefile
 *  DUMP_STRING_TO_FILE(..)
 *
 *  RESET_TIMER()     set timer to zero
 *  DUMP_TIMER()      dump time that has passed since the last
 *                    RESET_TIMER()
 */
#define DUMP_TIMER()    fprintf( DUMPOUT, "%s:%d Timer: %f\n", __FILE__, __LINE__, (float)(clock()-last_timer_reset)/(float)CLOCKS_PER_SEC )

/**
 *                         useful macros
 * -----------------------------------------------------------------------
 *  The following macros can be helpful to display variables during the
 *  runtime of the program.
 *  If you start the program in a shell from xemacs, you can jump to the
 *  respectove rows by switching to the compilation-minor-mode.
 *
 *  DUMP_POSITION()           dumps the actual position within the program
 *  DUMP_MESSAGE( .)          dump a message in addition
 *  DUMP_INT(..)              dump an integer variable
 *
 *  DUMP_MATRIX_TO_FILE(..)
 *  DUMP_INT_TO_FILE(..)      writes the value of the variable in
 *  DUMP_DOUBLE_TO_FILE(..)   a tracefile
 *  DUMP_STRING_TO_FILE(..)
 *
 *  RESET_TIMER()     set timer to zero
 *  DUMP_TIMER()      dump time that has passed since the last
 *                    RESET_TIMER()
 */
#define DUMP_MATRIX_TO_FILE( m, nrl, nrh, ncl, nch, xlength, ylength) \
        {  \
           static nCount = 0; \
	   char szFileName[100];  \
	   sprintf( szFileName, "%s__%d__%s.out", __FILE__, __LINE__, #m); \
           write_matrix( szFileName, m, nrl, nrh, ncl, nch, xlength, ylength, nCount == 0); \
	   ++nCount; \
        }

/**
 *                         useful macros
 * -----------------------------------------------------------------------
 *  The following macros can be helpful to display variables during the
 *  runtime of the program.
 *  If you start the program in a shell from xemacs, you can jump to the
 *  respectove rows by switching to the compilation-minor-mode.
 *
 *  DUMP_POSITION()           dumps the actual position within the program
 *  DUMP_MESSAGE( .)          dump a message in addition
 *  DUMP_INT(..)              dump an integer variable
 *
 *  DUMP_MATRIX_TO_FILE(..)
 *  DUMP_INT_TO_FILE(..)      writes the value of the variable in
 *  DUMP_DOUBLE_TO_FILE(..)   a tracefile
 *  DUMP_STRING_TO_FILE(..)
 *
 *  RESET_TIMER()     set timer to zero
 *  DUMP_TIMER()      dump time that has passed since the last
 *                    RESET_TIMER()
 */
#define DUMP_INT_TO_FILE(n) \
        {  \
           static nCount = 0; \
           FILE *fh = 0; \
	   char szFileName[100];  \
	   sprintf( szFileName, "%s__%d__%s.out", __FILE__, __LINE__, #n); \
	   if( nCount == 0) \
              fh = fopen( szFileName, "w"); \
           else  \
              fh = fopen( szFileName, "a"); \
           if( fh )  \
              fprintf( fh, "%d:%d\n", nCount, n ); \
           else  \
              THROW_ERROR("Fehler beim Dumpen");  \
           fclose(fh);  \
	   ++nCount; \
        }

/**
 *                         useful macros
 * -----------------------------------------------------------------------
 *  The following macros can be helpful to display variables during the
 *  runtime of the program.
 *  If you start the program in a shell from xemacs, you can jump to the
 *  respectove rows by switching to the compilation-minor-mode.
 *
 *  DUMP_POSITION()           dumps the actual position within the program
 *  DUMP_MESSAGE( .)          dump a message in addition
 *  DUMP_INT(..)              dump an integer variable
 *
 *  DUMP_MATRIX_TO_FILE(..)
 *  DUMP_INT_TO_FILE(..)      writes the value of the variable in
 *  DUMP_DOUBLE_TO_FILE(..)   a tracefile
 *  DUMP_STRING_TO_FILE(..)
 *
 *  RESET_TIMER()     set timer to zero
 *  DUMP_TIMER()      dump time that has passed since the last
 *                    RESET_TIMER()
 */
#define DUMP_DOUBLE_TO_FILE(d) \
        {  \
           static nCount = 0; \
           FILE *fh = 0; \
	   char szFileName[100];  \
	   sprintf( szFileName, "%s__%d__%s.out", __FILE__, __LINE__, #d); \
	   if( nCount == 0) \
              fh = fopen( szFileName, "w"); \
           else  \
              fh = fopen( szFileName, "a"); \
           if( fh )  \
              fprintf( fh, "%d:%f\n", nCount, d ); \
           else  \
              THROW_ERROR("Fehler beim Dumpen");  \
           fclose(fh);  \
	   ++nCount; \
        }

/**
 *                         useful macros
 * -----------------------------------------------------------------------
 *  The following macros can be helpful to display variables during the
 *  runtime of the program.
 *  If you start the program in a shell from xemacs, you can jump to the
 *  respectove rows by switching to the compilation-minor-mode.
 *
 *  DUMP_POSITION()           dumps the actual position within the program
 *  DUMP_MESSAGE( .)          dump a message in addition
 *  DUMP_INT(..)              dump an integer variable
 *
 *  DUMP_MATRIX_TO_FILE(..)
 *  DUMP_INT_TO_FILE(..)      writes the value of the variable in
 *  DUMP_DOUBLE_TO_FILE(..)   a tracefile
 *  DUMP_STRING_TO_FILE(..)
 *
 *  RESET_TIMER()     set timer to zero
 *  DUMP_TIMER()      dump time that has passed since the last
 *                    RESET_TIMER()
 */
#define DUMP_STRING_TO_FILE(s) \
        {  \
           static nCount = 0; \
           FILE *fh = 0; \
	   char szFileName[100];  \
	   sprintf( szFileName, "%s__%d__%s.out", __FILE__, __LINE__, #s); \
	   if( nCount == 0) \
              fh = fopen( szFileName, "w"); \
           else  \
              fh = fopen( szFileName, "a"); \
           if( fh )  \
              fprintf( fh, "%d:%s\n", nCount, s ); \
           else  \
              THROW_ERROR("Fehler beim Dumpen");  \
           fclose(fh);  \
	   ++nCount; \
        }

#endif     

