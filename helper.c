#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <sys/stat.h>
#include "helper.h"
#include "logger.h"
#include "boundary_configurator.h"
#include "boundary_val.h"



/* ----------------------------------------------------------------------- */
/*                             auxiliary functions                         */
/* ----------------------------------------------------------------------- */
int min( int a, int b)           
{
    if( a < b ) return a;
    return b;
}

int max( int a, int b)
{
    if( a > b ) return a;
    return b;
}

double fmin( double a, double b)
{
    if( a < b ) return a;
    return b;
}

double fmax( double a, double b)
{
    if( a > b ) return a;
    return b;
}
/* ----------------------------------------------------------------------- */
/*                             custom auxiliary functions                  */
/* ----------------------------------------------------------------------- */

// Returns 1 (True) if the cell is an obstacle
int isObstacle(int flag){
    return (flag>>CENTER)&1;
}

// Returns 1 (True) if the cell is fluid
int isFluid(int flag){
    return !((flag>>CENTER)&1);
}

// Returns 1 (True) if the cell is coupling
int isCoupling(int flag){
    return ((flag>>CBIT)&1);
}

// Returns 1 (True) if the neighbouring cell in the indicated direction is an obstacle
int isNeighbourObstacle(int flag, Direction direction){
    return (flag>>direction)&1;
}

// Returns 1 (True) if the neighbouring cell in the indicated direction is fluid
int isNeighbourFluid(int flag, Direction direction){
    return !((flag>>direction)&1);
}

// Returns 1 (True) if the cell is present at a corner (bordering only 2 fluid cells)
int isCorner(int flag){
//    return ((flag&(1<<TOP))>>TOP)^((flag&(1<<BOT))>>BOT) && ((flag&(1<<LEFT))>>LEFT)^((flag&(1<<RIGHT))>>RIGHT);
    return ((flag>>TOP)&1)^((flag>>BOT)&1) && ((flag>>LEFT)&1)^((flag>>RIGHT)&1) && isObstacle(flag);
}


// Computes skip condition for u-boundary value determination (if top, right and bottom cells are obstacles)
int skipU(int flag){
//    return (flag&(1<<TOP)) && (flag&(1<<RIGHT)) && (flag&(1<<BOT));
    return ((flag>>TOP)&1) && ((flag>>RIGHT)&1) && ((flag>>BOT)&1);
}

// Computes skip condition for v-boundary value determination (if left, top and right cells are obstacles)
int skipV(int flag){
//    return (flag&(1<<LEFT)) && (flag&(1<<TOP)) && (flag&(1<<RIGHT));
    return ((flag>>LEFT)&1) && ((flag>>TOP)&1) && ((flag>>RIGHT)&1);
}

// Function that checks geometry for forbidden cases
void geometryCheck(int** Flag, int imax, int jmax){
    int isForbidden = 0;
    for (int j = jmax; j > 0; j--)
    {
        for (int i = 1; i < imax + 1; i++)
        {
            if( isFluid(Flag[i][j]) )
            {
                continue;
            }
            else if( ( isFluid(Flag[i][j + 1]) && isFluid(Flag[i][j - 1]) ) ||
                 ( isFluid(Flag[i - 1][j]) && isFluid(Flag[i + 1][j]) )  )
            {
                logMsg(ERROR, "Forbidden Geometry present at (%d,%d)", i, j);
                isForbidden++;
            }
        }
    }
    if(isForbidden == 0){
        logMsg(PRODUCTION, "Geometry has no forbidden configurations!");
    }
    else
    {
        logMsg(ERROR, "%d forbidden geometries found!",isForbidden);
        THROW_ERROR("Forbidden geometries!");
    }
}


/* ----------------------------------------------------------------------- */
/*                         local auxiliary functions                       */
/* ----------------------------------------------------------------------- */

clock_t last_timer_reset;

int min_int( const int n1, const int n2 )
{
    if( n1 < n2 ) return n1;
    return n2;
}



/* ----------------------------------------------------------------------- */
/*                             read datafile                               */
/* ----------------------------------------------------------------------- */

void errhandler( int nLine, const char *szFile, const char *szString )
{
    int err = errno;

    fprintf( ERROUT, "%s:%d Error : %s", szFile, nLine, szString );
    fprintf( ERROUT, "\n" );
    
    /* if an error within the c-library occured, an error code can be   */
    /* found in the global variable err                                 */
    if( err != 0 )
    {
	fprintf( ERROUT, "C-Lib   errno    = %d\n", err);
	fprintf( ERROUT, "C-Lib   strerror = %s\n", strerror( err ) );
    }
    exit(1);
}


/*  for comfort */
#define READ_ERROR(szMessage, szVarName, szFileName, nLine) \
  { char szTmp[80]; \
    if( nLine ) \
	sprintf( szTmp, " %s  File: %s   Variable: %s  Line: %d", szMessage, szFileName, szVarName, nLine ); \
    else \
	sprintf( szTmp, " %s  File: %s   Variable: %s ", szMessage, szFileName, szVarName); \
    THROW_ERROR( szTmp ); \
  }
    

/* --------------------------------------------------------------------------*/
/* The function searches the datafile fh for the line defining the variable  */
/* szVarName and returns the respctive string including the value of the     */
/* variable. If there's no appropriate line within the datafile, the program */
/* stops with an error messsage.                                             */
/* ATTENTION: The pointer returned refers to a static variable within the    */
/* function. To maintain the string over several program calls, it has to be */
/* copied!!!                                                                 */
/*                                                                           */
char *find_string(const char *szFileName, const char *szVarName, Optional optional)
{ 
    int nLine = 0;
    int i;
    FILE *fh = NULL;
    
    static char szBuffer[MAX_LINE_LENGTH];	/* containes the line read  */
                                               /* from the datafile        */

    char* szLine = szBuffer;
    char* szValue = NULL;
    char* szName = NULL;

    /* open file */
    fh = fopen( szFileName, "rt" );
    if( fh == 0 ) 
	READ_ERROR("Could not open file", szVarName, szFileName, 0);

    /* searching */
    while( ! feof(fh) )
    {
	fgets( szLine, MAX_LINE_LENGTH, fh );
	++nLine;

	/* remove comments */
	for( i = 0; i < strlen(szLine); i++)
	    if( szLine[i] == '#' )
	    {
		szLine[i] = '\0'; /* Stringende setzen */
		break;
	    }

	/* remove empty lines */
	while( isspace( (int)*szLine ) && *szLine) ++szLine;
	if( strlen( szLine ) == 0) continue; 

	/* now, the name can be extracted */
	szName = szLine;
	szValue = szLine;
	while( (isalnum( (int)*szValue ) || *szValue == '_') && *szValue) ++szValue;
	
	/* is the value for the respective name missing? */
	if( *szValue == '\n' || strlen( szValue) == 0)  
	    READ_ERROR("wrong format", szName, szFileName, nLine);
	
	*szValue = 0;		/* complete szName! at the right place */
	++szValue;
        
	/* read next line if the correct name wasn't found */
	if( strcmp( szVarName, szName)) continue;

	/* remove all leading blnkets and tabs from the value string  */
	while( isspace( (int)*szValue) ) ++szValue;
	if( *szValue == '\n' || strlen( szValue) == 0)  
	    READ_ERROR("wrong format", szName, szFileName, nLine);
	
	fclose(fh);
	return szValue;
    }  
   
    if (optional == REQUIRED)
        READ_ERROR("variable not found", szVarName, szFileName, nLine);
    
    return NULL;		/* dummy to satisfy the compiler  */
} 

void read_string(const char *szFileName, const char *szVarName, char *pVariable, Optional optional)
{
    char* szValue = NULL;	/* string containg the read variable value */

    if( szVarName  == 0 )  THROW_ERROR("null pointer given as variable name" );
    if( szFileName == 0 )  THROW_ERROR("null pointer given as filename" );
    if( pVariable  == 0 )  THROW_ERROR("null pointer given as variable" );

    if( szVarName[0] == '*' )
	szValue = find_string(szFileName, szVarName + 1, optional);
    else
	szValue = find_string(szFileName, szVarName, optional);
    
    if (szValue)
    {
        if (sscanf(szValue, "%s", pVariable) == 0)
        READ_ERROR("wrong format", szVarName, szFileName, 0);
    }
    else // If not found default to 0
        strcpy(pVariable, "NULLSTRING");
    
    logMsg(PRODUCTION, "File: %s\t\t%s%s= %s", szFileName,
            szVarName,
            &("               "[min_int( strlen(szVarName), 15)]),
            pVariable );
}

void read_int(const char *szFileName, const char *szVarName, int *pVariable, Optional optional)
{
    char* szValue = NULL;	/* string containing the read variable value */

    if( szVarName  == 0 )  THROW_ERROR("null pointer given as varable name" );
    if( szFileName == 0 )  THROW_ERROR("null pointer given as filename" );
    if( pVariable  == 0 )  THROW_ERROR("null pointer given as variable" );

    if( szVarName[0] == '*' )
	szValue = find_string(szFileName, szVarName + 1, optional);
    else
	szValue = find_string(szFileName, szVarName, optional);
    
    if (szValue)
    {
        if (sscanf(szValue, "%d", pVariable) == 0)
        READ_ERROR("wrong format", szVarName, szFileName, 0);
    }
    else // If not found default to 0
        *pVariable = 0;
    
    logMsg(PRODUCTION, "File: %s\t\t%s%s= %d", szFileName,
            szVarName,
            &("               "[min_int( strlen(szVarName), 15)]),
            *pVariable );
}

void read_double(const char *szFileName, const char *szVarName, double *pVariable, Optional optional)
{
    char* szValue = NULL;	/* String mit dem eingelesenen Variablenwert */

    if( szVarName  == 0 )  THROW_ERROR("null pointer given as varable name" );
    if( szFileName == 0 )  THROW_ERROR("null pointer given as filename" );
    if( pVariable  == 0 )  THROW_ERROR("null pointer given as variable" );

    if( szVarName[0] == '*' )
	szValue = find_string(szFileName, szVarName + 1, optional);
    else
        szValue = find_string(szFileName, szVarName, optional);
    
    if (szValue)
    {
        if (sscanf(szValue, "%lf", pVariable) == 0)
        READ_ERROR("wrong format", szVarName, szFileName, 0);
    }
    else // If not found default to 0
        *pVariable = 0.0;
    
    logMsg(PRODUCTION, "File: %s\t\t%s%s= %f", szFileName,
            szVarName,
            &("               "[min_int( strlen(szVarName), 15)]),
            *pVariable );
}


/* ----------------------------------------------------------------------- */
/*                   write matrices to a file                              */
/* ----------------------------------------------------------------------- */

void write_matrix( const char* szFileName,       /* filename */
		   double **m,		       /* matrix */
		   int nrl,		       /* first column */
		   int nrh,		       /* last column */
		   int ncl,		       /* first row */
		   int nch,		       /* last row */
		 double xlength,	       /* size of the geometry in */
                                               /* x-direction */
		 double ylength,	       /* size of the geometry in */
                                               /* y-direction  */
		   int fFirst ) 	       /* 0 == append, else overwrite*/
{
   int i, j;
   FILE * fh = 0;
   int nSize = (nrh-nrl+1) * (nch-ncl+1);
   float *tmp = (float *)malloc( (size_t)(nSize * sizeof(float)));
   int k = 0;

   if( fFirst )				/* first call of the function ? */
   {
       fh = fopen( szFileName, "w");	/* overwrite file/write new file */
       if( fh == NULL )			/* opening failed ? */
       {
	   char szBuff[80];
	   sprintf( szBuff, "Outputfile %s cannot be created", szFileName );
	   THROW_ERROR( szBuff );
       }
       
/*       fprintf( fh,"%f\n%f\n%d\n%d\n%d\n%d\n", xlength, ylength, nrl, nrh, ncl, nch ); */
   }
   else
   {
       fh = fopen( szFileName ,"a");	/* append to the file */
       if( fh == NULL )			/* opening failed ? */
       {
	   char szBuff[80];
	   sprintf( szBuff, "Outputfile %s cannot be opened", szFileName );
	   THROW_ERROR( szBuff );
       }
   } 

   for( j = ncl; j <= nch; j++)
       for( i = nrl; i <= nrh; i++)
	   tmp[k++] = (float)m[i][j];

   fwrite( tmp, sizeof(float), nSize, fh);

   if( fclose(fh) )
   {
       char szBuff[80];
       sprintf( szBuff, "Outputfile %s cannot be closed", szFileName );
       THROW_ERROR( szBuff );
   };

   free( tmp );
}


void read_matrix( const char* szFileName,       /* filename */
		   double **m,		       /* matrix */
		   int nrl,		       /* first column */
		   int nrh,		       /* last column */
		   int ncl,		       /* first row */
		   int nch		       /* last row */
                  ) 	  
{
   int i, j;
   FILE * fh = 0;
   int nSize = (nrh-nrl+1) * (nch-ncl+1);
   float *tmp = (float *)malloc( (size_t)(nSize * sizeof(float)));
   int k = 0;

       fh = fopen( szFileName, "r");	/* overwrite file/write new file */
       if( fh == NULL )			/* opening failed ? */
       {
	   char szBuff[80];
	   sprintf( szBuff, "Can not read file %s !!!", szFileName );
	   THROW_ERROR( szBuff );
       }


   fread( tmp, sizeof(float), nSize, fh);

   for( j = ncl; j <= nch; j++)
       for( i = nrl; i <= nrh; i++)
	   m[i][j]=tmp[k++];

   if( fclose(fh) )
   {
       char szBuff[80];
       /*orig bug:
       sscanf( szBuff, "Inputfile %s cannot be closed", szFileName );*/
       sprintf( szBuff, "Inputfile %s cannot be closed", szFileName );
       THROW_ERROR( szBuff );
   };

   free( tmp );
}


/* ----------------------------------------------------------------------- */
/*                      general matrix functions                           */
/* ----------------------------------------------------------------------- */

/*  allocates storage for a matrix                                         */
double **matrix( int nrl, int nrh, int ncl, int nch )
{
   int i;
   int nrow = nrh - nrl + 1;	/* compute number of lines */
   int ncol = nch - ncl + 1;	/* compute number of columns */
   
   double **pArray  = (double **) malloc((size_t)( nrow * sizeof(double*)) );
   double  *pMatrix = (double *)  malloc((size_t)( nrow * ncol * sizeof( double )));

   if( pArray  == 0)  THROW_ERROR("Storage cannot be allocated");
   if( pMatrix == 0)  THROW_ERROR("Storage cannot be allocated");

   /* first entry of the array points to the value corrected by the 
      beginning of the column */
   pArray[0] = pMatrix - ncl; 

   /* compute the remaining array entries */
   for( i = 1; i < nrow; i++ )
   {
       pArray[i] = pArray[i-1] + ncol;
   }

   /* return the value corrected by the beginning of a line */
   return pArray - nrl;
}


/* deallocates the storage of a matrix  */
void free_matrix( double **m, int nrl, int nrh, int ncl, int nch )
{
   double **pArray  = m + nrl;
   double  *pMatrix = m[nrl]+ncl;

   free( pMatrix );
   free( pArray );
}

void init_matrix( double **m, int nrl, int nrh, int ncl, int nch, double a)
{
   int i,j;
   for( i = nrl; i <= nrh; i++)
       for( j = ncl; j <= nch; j++)
	   m[i][j] = a;
}


/* allocates storage for a matrix */
int **imatrix( int nrl, int nrh, int ncl, int nch )
{
   int i;

   int nrow = nrh - nrl + 1;	/* compute number of rows */
   int ncol = nch - ncl + 1;	/* compute number of columns */
   
   int **pArray  = (int **) malloc((size_t)( nrow * sizeof( int* )) );
   int  *pMatrix = (int *)  malloc((size_t)( nrow * ncol * sizeof( int )));


   if( pArray  == 0)  THROW_ERROR("Storage cannot be allocated");
   if( pMatrix == 0)  THROW_ERROR("Storage cannot be allocated");

   /* first entry of the array points to the value corrected by the 
      beginning of the column */
   pArray[0] = pMatrix - ncl; 

   /* compute the remaining array entries */
   for( i = 1; i < nrow; i++ )
   {
       pArray[i] = pArray[i-1] + ncol;
   }

   /* return the value corrected by the beginning of a line */
   return pArray - nrl;
}

/* deallocates the storage of a matrix  */
void free_imatrix( int **m, int nrl, int nrh, int ncl, int nch )
{
   int **pArray  = m + nrl;
   int  *pMatrix = m[nrl]+ncl;

   free( pMatrix );
   free( pArray );
}

void init_imatrix( int **m, int nrl, int nrh, int ncl, int nch, int a)
{
   int i,j;
   for( i = nrl; i <= nrh; i++)
       for( j = ncl; j <= nch; j++)
	   m[i][j] = a;
}


int **read_pgm(const char *filename) {
    logMsg(DEBUG, "read_pgm(): filename = %s", filename); //debug
    FILE *input = NULL;
    char line[1024];
    int levels;
    int xsize, ysize;
    int **pic = NULL;

    if ((input = fopen(filename, "rb")) == 0) {
        char szBuff[80];
        sprintf(szBuff, "Can not read file %s !!!", filename);
        THROW_ERROR(szBuff);
    }

    /* check for the right "magic number" */
    if (fread(line, 1, 3, input) != 3) {
        fclose(input);
        THROW_ERROR("Error Wrong Magic field!");
    }

    /* skip the comments */
    do
        fgets(line, sizeof line, input);
    while (*line == '#');

    /* read the width and height */
    sscanf(line, "%d %d\n", &xsize, &ysize);

    logMsg(PRODUCTION, "Image size: %d x %d", xsize, ysize);

    /* read # of gray levels */
    fgets(line, sizeof line, input);
    sscanf(line, "%d\n", &levels);

    /* allocate memory for image */
    pic = imatrix(0, xsize-1, 0, ysize-1);
    logMsg(PRODUCTION, "Image initialised...");

    for (int j = ysize-1; j >= 0; --j) {
        for (int i = 0; i < xsize; ++i) {
            int value;
            fscanf(input, "%d", &value);

            if (value == EOF) {
                fclose(input);
                THROW_ERROR("read of geometry file failed!");
            }
            pic[i][j] = value;
        }
    }

    /* close file */
    fclose(input);

    return pic;
}

/* decodes Flag back into *.pgm format */
void decode_flags(int imax, int jmax, int **Flag, int** pic)
{
    // always running in EXTENDED mode


        // decode flags into gray colour map - go through the whole domain including boundaries, i.e. EXTENDED mode
        for (int i = 0; i <= imax+1; i++)
        {
            for (int j = 0; j <= jmax+1; j++)
            {
                pic[i][j] = (((1 << CENTER) & Flag[i][j]) == 0) * FLUID_PIXEL // if the cell is fluid then set pic value to FLUID_PIXEL
                            + (((1 << NSBIT) & Flag[i][j]) != 0) * NOSLIP_PIXEL // ...
                            + (((1 << FSBIT) & Flag[i][j]) != 0) * FREESLIP_PIXEL // ...
                            + (((1 << OFBIT) & Flag[i][j]) != 0) * OUTFLOW_PIXEL // ...
                            + (((1 << IFBIT) & Flag[i][j]) != 0) * INFLOW_PIXEL // ...
                            + (((1 << CBIT) & Flag[i][j]) != 0) * COUPLING_PIXEL; // if the cell is coupling then set pic value to COUPLING_PIXEXL
            }
        }
}


void write_pgm(int xsize, int ysize, int **pgm, const char *outputFolder, const char *szProblem, int iterationNumber)
{
    // always running in EXTENDED mode
    char szFileName[80];
    char szFileFullPath[512];
    FILE *fp = NULL;
    sprintf(szFileName, "%s.%i.pgm", szProblem, iterationNumber);
    sprintf(szFileFullPath, "%s/%s", outputFolder, szFileName);
    fp = fopen(szFileFullPath, "wb");
    if (fp == NULL)
    {
        char szBuff[80];
        logMsg(ERROR, "Failed to open %s", szFileFullPath);
        sprintf(szBuff, "Failed to open %s", szFileName);
        THROW_ERROR(szBuff);
        return;
    }


        fprintf(fp, "P2\n");
        //printf("P2\n");

        fprintf(fp, "# 0 - no-slip\n");
        //printf("# 0 - no-slip\n");

        fprintf(fp, "# 1 - free-slip\n");
        //printf("# 1 - free-slip\n");

        fprintf(fp, "# 2 - outflow\n");
        //printf("# 2 - outflow\n");

        fprintf(fp, "# 3 - inflow\n");
        //printf("# 3 - inflow\n");

        fprintf(fp, "# 4 - coupling\n");
        //printf("# 4 - coupling\n");

        fprintf(fp, "# 6 - fluid\n");
        //printf("# 6 - fluid\n");

        fprintf(fp, "%d %d\n", xsize, ysize);
        //printf("%d %d\n", xsize, ysize);

        fprintf(fp, "%d\n", FLUID_PIXEL);
        //printf("%d\n", FLUID_PIXEL);

        // create buffer - one row of *.pgm picture
        char* buffer = (char*)malloc((xsize*2+1) * sizeof(char));

        for (int j = ysize -1; j >= 0; j--)
        {
            int i = 0;

            for(i = 0; i <= 2*xsize-1; i = i+2)
            {
                // fill in line buffer including spaces inbetween the numbers
                buffer[i] = 48 + pgm[i/2][j];
                buffer[i + 1] = ' ';
            }

            // add newline character
            buffer[2*xsize-1] = '\n';
            // add end of string character
            buffer[2*xsize] = '\0';

            fprintf(fp, "%s", buffer);
            //printf("%s\n",buffer);
        }

        // on Windows sth strange here - infinite loop (?) how to correctly free the buffer (?)
        free(buffer);
        buffer = NULL;
        fclose(fp);
        //printf("EOF");
}


void flipToFluid(double **U, double **V, int  **Flag, int i, int j)
{
    if (isObstacle(Flag[i][j])){ printf("eeeerrrrrooooorrr\n");}

    Flag[i][j] = 0
                  + (1 << TOP) * isObstacle(Flag[i][j + 1])
                  + (1 << BOT) * isObstacle(Flag[i][j - 1])
                  + (1 << LEFT) * isObstacle(Flag[i - 1][j])
                  + (1 << RIGHT) * isObstacle(Flag[i + 1][j])
                 + (1 << TBIT) * (NEUMANN); // NEUMANN by default

    Flag[i+1][j] -= (1 << LEFT);
    Flag[i-1][j] -= (1 << RIGHT);
    Flag[i][j-1] -= (1 << TOP);
    Flag[i][j+1] -= (1 << BOT);

    U[i][j] = 0;
    V[i][j] = 0;

}

void flipToSolid(double **U, double **V, double** P, int  **Flag, int i, int j)
{
    Flag[i][j] = (1 << CENTER)
                 + (1 << NSBIT)
                 + (1 << TOP) * isObstacle(Flag[i][j + 1])
                 + (1 << BOT) * isObstacle(Flag[i][j - 1])
                 + (1 << LEFT) * isObstacle(Flag[i - 1][j])
                 + (1 << RIGHT) * isObstacle(Flag[i + 1][j])
    + (1 << TBIT) * (NEUMANN); // NEUMANN by default

    Flag[i+1][j] += (1 << LEFT);
    Flag[i-1][j] += (1 << RIGHT);
    Flag[i][j -1] += (1 << TOP);
    Flag[i][j + 1] += (1 << BOT);

    U[i][j] = 0;
    V[i][j] = 0;
    P[i][j] = 0;
}

int checkVelocityMagnitude(double eps, double U, double V)
{
    if(sqrt( pow(U,2) + pow(V,2) ) < eps)
        return 1;
    else
        return 0;

}

void update_pgm(int imax, int jmax, int *noFluidCells, int **pgm, double **Flag, double **P, double **U, double **V, double eps)
{
    //int i = 1;
    //int j = 1;
    int isFlip = 0;

    for (int i = 1; i < imax + 1; i++)
    {
        for(int j = 1; j < jmax + 1; j++)
        {
            isFlip = 0;
            if (isObstacle(Flag[i][j]))
            {
                //isFlip = checkVelocityMagnitude(eps,U[i][j],V[i][j]);
                if (isFlip) {
                    printf("inside flip solid\n");
                    flipToFluid(U,V, Flag, i, j);
                    (*noFluidCells)++;
                }
            }
            else
            {
                isFlip = checkVelocityMagnitude(eps,U[i][j],V[i][j]);
                if (isFlip)
                {
                    flipToSolid(U,V, P, Flag, i, j);
                    (*noFluidCells)--;
                }
            }

        }
    }

}