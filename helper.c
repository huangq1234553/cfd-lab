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
int min(int a, int b)
{
    if (a < b) return a;
    return b;
}

int max(int a, int b)
{
    if (a > b) return a;
    return b;
}

double fmin(double a, double b)
{
    if (a < b) return a;
    return b;
}

double fmax(double a, double b)
{
    if (a > b) return a;
    return b;
}

/* ----------------------------------------------------------------------- */
/*                             custom auxiliary functions                  */
/* ----------------------------------------------------------------------- */

// Returns 1 (True) if the cell is an obstacle
inline int isObstacle(int flag)
{
    return (flag >> CENTER) & 1;
}

// Returns 1 (True) if the cell is fluid
inline int isFluid(int flag)
{
    return !((flag >> CENTER) & 1);
}

// Returns 1 (True) if the cell is coupling
inline int isCoupling(int flag)
{
    return ((flag >> CBIT) & 1);
}

// Returns 1 (True) if the cell is outflow
inline int isOutflow(int flag)
{
    return ((flag >> OFBIT) & 1);
}

// Returns 1 (True) if the cell is inflow
inline int isInflow(int flag)
{
    return ((flag >> IFBIT) & 1);
}

// Returns 1 (True) if cell cannot be switched between solid & fluid
inline int isGeometryConstant(int flag)
{
    return ((flag >> GEOMETRYMASKBIT) & 1);
}

// Returns 1 (True) if the neighbouring cell in the indicated direction is an obstacle
inline int isNeighbourObstacle(int flag, Direction direction)
{
    return (flag >> direction) & 1;
}

// Returns 1 (True) if the neighbouring cell in the indicated direction is fluid
inline int isNeighbourFluid(int flag, Direction direction)
{
    return !((flag >> direction) & 1);
}

// Current cell's neighbor in the specified direction allows current to be flipped to solid
int doesNeighbourAllowFlipToSolid(int **Flags, int i, int j, Direction direction, int imax, int jmax)
{
    int cell = Flags[i][j];
    int neighbour = 0;
    int iN, jN;
    switch (direction)
    {
        case LEFT:
            iN = i - 1;
            jN = j;
            break;
        case RIGHT:
            iN = i + 1;
            jN = j;
            break;
        case TOP:
            iN = i;
            jN = j + 1;
            break;
        case BOT:
            iN = i;
            jN = j - 1;
            break;
        default:
            iN = i;
            jN = j;
    }
    neighbour = Flags[iN][jN];
    int check = (iN > 0) && (iN <= imax) && (jN > 0) && (jN <= jmax);
    return check && isNeighbourObstacle(cell, direction) && !isOutflow(neighbour) && !isInflow(neighbour);
}

inline int hasAtLeastOneObstacleNeighbour(int flag)
{
    return isNeighbourObstacle(flag, LEFT)
           || isNeighbourObstacle(flag, RIGHT)
           || isNeighbourObstacle(flag, BOT)
           || isNeighbourObstacle(flag, TOP);
}

// Returns 1 (True) if the cell is present at a corner (bordering only 2 fluid cells)
inline int isCorner(int flag)
{
//    return ((flag&(1<<TOP))>>TOP)^((flag&(1<<BOT))>>BOT) && ((flag&(1<<LEFT))>>LEFT)^((flag&(1<<RIGHT))>>RIGHT);
    return ((flag >> TOP) & 1) ^ ((flag >> BOT) & 1) && ((flag >> LEFT) & 1) ^ ((flag >> RIGHT) & 1) &&
           isObstacle(flag);
}


// Computes skip condition for u-boundary value determination (if top, right and bottom cells are obstacles)
inline int skipU(int flag)
{
//    return (flag&(1<<TOP)) && (flag&(1<<RIGHT)) && (flag&(1<<BOT));
    return ((flag >> TOP) & 1) && ((flag >> RIGHT) & 1) && ((flag >> BOT) & 1);
}

// Computes skip condition for v-boundary value determination (if left, top and right cells are obstacles)
inline int skipV(int flag)
{
//    return (flag&(1<<LEFT)) && (flag&(1<<TOP)) && (flag&(1<<RIGHT));
    return ((flag >> LEFT) & 1) && ((flag >> TOP) & 1) && ((flag >> RIGHT) & 1);
}

// Function that checks geometry for forbidden cases
void geometryCheck(int **Flags, int imax, int jmax, int *noFluidCells, bool fixInitialGeometry)
{
    int isForbidden = 0;
    for (int j = jmax; j > 0; j--)
    {
        for (int i = 1; i < imax + 1; i++)
        {
            if (isFluid(Flags[i][j]))
            {
                continue;
            }
            else if ((isFluid(Flags[i][j + 1]) && isFluid(Flags[i][j - 1])) ||
                     (isFluid(Flags[i - 1][j]) && isFluid(Flags[i + 1][j])))
            {
                logMsg(ERROR, "Forbidden Geometry present at (%d,%d)", i, j);
                isForbidden++;
            }
        }
    }
    if (isForbidden == 0)
    {
        logMsg(PRODUCTION, "Geometry has no forbidden configurations!");
    }
    else
    {
        logMsg(ERROR, "%d forbidden geometries found!", isForbidden);
        if (fixInitialGeometry)
        {
            geometryFix(NULL, NULL, NULL, Flags, imax, jmax, noFluidCells, NULL);
            logMsg(PRODUCTION, "Geometry fixed! New total fluid cells in domain: %d", (*noFluidCells));
        }
        else
        {
            THROW_ERROR("Forbidden geometries!");
        }
    }
}


/* ----------------------------------------------------------------------- */
/*                         local auxiliary functions                       */
/* ----------------------------------------------------------------------- */

clock_t last_timer_reset;

int min_int(const int n1, const int n2)
{
    if (n1 < n2) return n1;
    return n2;
}



/* ----------------------------------------------------------------------- */
/*                             read datafile                               */
/* ----------------------------------------------------------------------- */

void errhandler(int nLine, const char *szFile, const char *szString)
{
    int err = errno;
    
    fprintf(ERROUT, "%s:%d Error : %s", szFile, nLine, szString);
    fprintf(ERROUT, "\n");
    
    /* if an error within the c-library occured, an error code can be   */
    /* found in the global variable err                                 */
    if (err != 0)
    {
        fprintf(ERROUT, "C-Lib   errno    = %d\n", err);
        fprintf(ERROUT, "C-Lib   strerror = %s\n", strerror(err));
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
    
    static char szBuffer[MAX_LINE_LENGTH];    /* containes the line read  */
    /* from the datafile        */
    
    char *szLine = szBuffer;
    char *szValue = NULL;
    char *szName = NULL;
    
    /* open file */
    fh = fopen(szFileName, "rt");
    if (fh == 0)
    READ_ERROR("Could not open file", szVarName, szFileName, 0);
    
    /* searching */
    while (!feof(fh))
    {
        fgets(szLine, MAX_LINE_LENGTH, fh);
        ++nLine;
        
        /* remove comments */
        for (i = 0; i < strlen(szLine); i++)
        {
            if (szLine[i] == '#')
            {
                szLine[i] = '\0'; /* Stringende setzen */
                break;
            }
        }
        
        /* remove empty lines */
        while (isspace((int) *szLine) && *szLine) ++szLine;
        if (strlen(szLine) == 0) continue;
        
        /* now, the name can be extracted */
        szName = szLine;
        szValue = szLine;
        while ((isalnum((int) *szValue) || *szValue == '_') && *szValue) ++szValue;
        
        /* is the value for the respective name missing? */
        if (*szValue == '\n' || strlen(szValue) == 0)
        READ_ERROR("wrong format", szName, szFileName, nLine);
        
        *szValue = 0;        /* complete szName! at the right place */
        ++szValue;
        
        /* read next line if the correct name wasn't found */
        if (strcmp(szVarName, szName)) continue;
        
        /* remove all leading blnkets and tabs from the value string  */
        while (isspace((int) *szValue)) ++szValue;
        if (*szValue == '\n' || strlen(szValue) == 0)
        READ_ERROR("wrong format", szName, szFileName, nLine);
        
        fclose(fh);
        return szValue;
    }
    
    if (optional == REQUIRED)
    READ_ERROR("variable not found", szVarName, szFileName, nLine);
    
    return NULL;        /* dummy to satisfy the compiler  */
}

void read_string(const char *szFileName, const char *szVarName, char *pVariable, Optional optional)
{
    char *szValue = NULL;    /* string containg the read variable value */
    
    if (szVarName == 0) THROW_ERROR("null pointer given as variable name");
    if (szFileName == 0) THROW_ERROR("null pointer given as filename");
    if (pVariable == 0) THROW_ERROR("null pointer given as variable");
    
    if (szVarName[0] == '*')
    {
        szValue = find_string(szFileName, szVarName + 1, optional);
    }
    else
    {
        szValue = find_string(szFileName, szVarName, optional);
    }
    
    if (szValue)
    {
        if (sscanf(szValue, "%s", pVariable) == 0)
        READ_ERROR("wrong format", szVarName, szFileName, 0);
    }
    else
    { // If not found default to 0
        strcpy(pVariable, "NULLSTRING");
    }
    
    logMsg(PRODUCTION, "File: %s\t\t%s%s= %s", szFileName,
           szVarName,
           &("               "[min_int(strlen(szVarName), 15)]),
           pVariable);
}

void read_int(const char *szFileName, const char *szVarName, int *pVariable, Optional optional)
{
    char *szValue = NULL;    /* string containing the read variable value */
    
    if (szVarName == 0) THROW_ERROR("null pointer given as varable name");
    if (szFileName == 0) THROW_ERROR("null pointer given as filename");
    if (pVariable == 0) THROW_ERROR("null pointer given as variable");
    
    if (szVarName[0] == '*')
    {
        szValue = find_string(szFileName, szVarName + 1, optional);
    }
    else
    {
        szValue = find_string(szFileName, szVarName, optional);
    }
    
    if (szValue)
    {
        if (sscanf(szValue, "%d", pVariable) == 0)
        READ_ERROR("wrong format", szVarName, szFileName, 0);
    }
    else
    { // If not found default to 0
        *pVariable = 0;
    }
    
    logMsg(PRODUCTION, "File: %s\t\t%s%s= %d", szFileName,
           szVarName,
           &("               "[min_int(strlen(szVarName), 15)]),
           *pVariable);
}

void read_double(const char *szFileName, const char *szVarName, double *pVariable, Optional optional)
{
    char *szValue = NULL;    /* String mit dem eingelesenen Variablenwert */
    
    if (szVarName == 0) THROW_ERROR("null pointer given as varable name");
    if (szFileName == 0) THROW_ERROR("null pointer given as filename");
    if (pVariable == 0) THROW_ERROR("null pointer given as variable");
    
    if (szVarName[0] == '*')
    {
        szValue = find_string(szFileName, szVarName + 1, optional);
    }
    else
    {
        szValue = find_string(szFileName, szVarName, optional);
    }
    
    if (szValue)
    {
        if (sscanf(szValue, "%lf", pVariable) == 0)
        READ_ERROR("wrong format", szVarName, szFileName, 0);
    }
    else
    { // If not found default to 0
        *pVariable = 0.0;
    }
    
    logMsg(PRODUCTION, "File: %s\t\t%s%s= %f", szFileName,
           szVarName,
           &("               "[min_int(strlen(szVarName), 15)]),
           *pVariable);
}


/* ----------------------------------------------------------------------- */
/*                   write matrices to a file                              */
/* ----------------------------------------------------------------------- */

void write_matrix(const char *szFileName,       /* filename */
                  double **m,               /* matrix */
                  int nrl,               /* first column */
                  int nrh,               /* last column */
                  int ncl,               /* first row */
                  int nch,               /* last row */
                  double xlength,           /* size of the geometry in */
        /* x-direction */
                  double ylength,           /* size of the geometry in */
        /* y-direction  */
                  int fFirst)           /* 0 == append, else overwrite*/
{
    int i, j;
    FILE *fh = 0;
    int nSize = (nrh - nrl + 1) * (nch - ncl + 1);
    float *tmp = (float *) malloc((size_t) (nSize * sizeof(float)));
    int k = 0;
    
    if (fFirst)                /* first call of the function ? */
    {
        fh = fopen(szFileName, "w");    /* overwrite file/write new file */
        if (fh == NULL)            /* opening failed ? */
        {
            char szBuff[80];
            sprintf(szBuff, "Outputfile %s cannot be created", szFileName);
            THROW_ERROR(szBuff);
        }

/*       fprintf( fh,"%f\n%f\n%d\n%d\n%d\n%d\n", xlength, ylength, nrl, nrh, ncl, nch ); */
    }
    else
    {
        fh = fopen(szFileName, "a");    /* append to the file */
        if (fh == NULL)            /* opening failed ? */
        {
            char szBuff[80];
            sprintf(szBuff, "Outputfile %s cannot be opened", szFileName);
            THROW_ERROR(szBuff);
        }
    }
    
    for (j = ncl; j <= nch; j++)
    {
        for (i = nrl; i <= nrh; i++)
        {
            tmp[k++] = (float) m[i][j];
        }
    }
    
    fwrite(tmp, sizeof(float), nSize, fh);
    
    if (fclose(fh))
    {
        char szBuff[80];
        sprintf(szBuff, "Outputfile %s cannot be closed", szFileName);
        THROW_ERROR(szBuff);
    };
    
    free(tmp);
}


void read_matrix(const char *szFileName,       /* filename */
                 double **m,               /* matrix */
                 int nrl,               /* first column */
                 int nrh,               /* last column */
                 int ncl,               /* first row */
                 int nch               /* last row */
)
{
    int i, j;
    FILE *fh = 0;
    int nSize = (nrh - nrl + 1) * (nch - ncl + 1);
    float *tmp = (float *) malloc((size_t) (nSize * sizeof(float)));
    int k = 0;
    
    fh = fopen(szFileName, "r");    /* overwrite file/write new file */
    if (fh == NULL)            /* opening failed ? */
    {
        char szBuff[80];
        sprintf(szBuff, "Can not read file %s !!!", szFileName);
        THROW_ERROR(szBuff);
    }
    
    
    fread(tmp, sizeof(float), nSize, fh);
    
    for (j = ncl; j <= nch; j++)
    {
        for (i = nrl; i <= nrh; i++)
        {
            m[i][j] = tmp[k++];
        }
    }
    
    if (fclose(fh))
    {
        char szBuff[80];
        /*orig bug:
        sscanf( szBuff, "Inputfile %s cannot be closed", szFileName );*/
        sprintf(szBuff, "Inputfile %s cannot be closed", szFileName);
        THROW_ERROR(szBuff);
    };
    
    free(tmp);
}


/* ----------------------------------------------------------------------- */
/*                      general matrix functions                           */
/* ----------------------------------------------------------------------- */

/*  allocates storage for a matrix                                         */
double **matrix(int nrl, int nrh, int ncl, int nch)
{
    int i;
    int nrow = nrh - nrl + 1;    /* compute number of lines */
    int ncol = nch - ncl + 1;    /* compute number of columns */
    
    double **pArray = (double **) malloc((size_t) (nrow * sizeof(double *)));
    double *pMatrix = (double *) malloc((size_t) (nrow * ncol * sizeof(double)));
    
    if (pArray == 0) THROW_ERROR("Storage cannot be allocated");
    if (pMatrix == 0) THROW_ERROR("Storage cannot be allocated");
    
    /* first entry of the array points to the value corrected by the
       beginning of the column */
    pArray[0] = pMatrix - ncl;
    
    /* compute the remaining array entries */
    for (i = 1; i < nrow; i++)
    {
        pArray[i] = pArray[i - 1] + ncol;
    }
    
    /* return the value corrected by the beginning of a line */
    return pArray - nrl;
}


/* deallocates the storage of a matrix  */
void free_matrix(double **m, int nrl, int nrh, int ncl, int nch)
{
    double **pArray = m + nrl;
    double *pMatrix = m[nrl] + ncl;
    
    free(pMatrix);
    free(pArray);
}

void init_matrix(double **m, int nrl, int nrh, int ncl, int nch, double a)
{
    int i, j;
    for (i = nrl; i <= nrh; i++)
    {
        for (j = ncl; j <= nch; j++)
        {
            m[i][j] = a;
        }
    }
}


/* allocates storage for a matrix */
int **imatrix(int nrl, int nrh, int ncl, int nch)
{
    int i;
    
    int nrow = nrh - nrl + 1;    /* compute number of rows */
    int ncol = nch - ncl + 1;    /* compute number of columns */
    
    int **pArray = (int **) calloc((size_t) nrow, sizeof(int *));
    int *pMatrix = (int *) calloc((size_t) nrow * ncol, sizeof(int));
    
    
    if (pArray == 0) THROW_ERROR("Storage cannot be allocated");
    if (pMatrix == 0) THROW_ERROR("Storage cannot be allocated");
    
    /* first entry of the array points to the value corrected by the
       beginning of the column */
    pArray[0] = pMatrix - ncl;
    
    /* compute the remaining array entries */
    for (i = 1; i < nrow; i++)
    {
        pArray[i] = pArray[i - 1] + ncol;
    }
    
    /* return the value corrected by the beginning of a line */
    return pArray - nrl;
}

/* deallocates the storage of a matrix  */
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch)
{
    int **pArray = m + nrl;
    int *pMatrix = m[nrl] + ncl;
    
    free(pMatrix);
    free(pArray);
}

void init_imatrix(int **m, int nrl, int nrh, int ncl, int nch, int a)
{
    int i, j;
    for (i = nrl; i <= nrh; i++)
    {
        for (j = ncl; j <= nch; j++)
        {
            m[i][j] = a;
        }
    }
}


int **read_pgm(const char *filename)
{
    logMsg(DEBUG, "read_pgm(): filename = %s", filename); //debug
    FILE *input = NULL;
    char line[1024];
    int levels;
    int xsize, ysize;
    int **pic = NULL;
    
    if ((input = fopen(filename, "rb")) == 0)
    {
        char szBuff[80];
        sprintf(szBuff, "Can not read file %s !!!", filename);
        THROW_ERROR(szBuff);
    }
    
    /* check for the right "magic number" */
    if (fread(line, 1, 3, input) != 3)
    {
        fclose(input);
        THROW_ERROR("Error Wrong Magic field!");
    }
    
    /* skip the comments */
    do
    {
        fgets(line, sizeof line, input);
    } while (*line == '#');
    
    /* read the width and height */
    sscanf(line, "%d %d\n", &xsize, &ysize);
    
    logMsg(PRODUCTION, "Image size: %d x %d", xsize, ysize);
    
    /* read # of gray levels */
    fgets(line, sizeof line, input);
    sscanf(line, "%d\n", &levels);
    
    /* allocate memory for image */
    pic = imatrix(0, xsize - 1, 0, ysize - 1);
    logMsg(PRODUCTION, "Image initialised...");
    
    for (int j = ysize - 1; j >= 0; --j)
    {
        for (int i = 0; i < xsize; ++i)
        {
            int value;
            fscanf(input, "%d", &value);
            
            if (value == EOF)
            {
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
void decode_flags(int imax, int jmax, int **Flag, int **pic)
{
    // always running in EXTENDED mode
    
    // decode flags into gray colour map - go through the whole domain including boundaries, i.e. EXTENDED mode
    for (int i = 0; i <= imax + 1; i++)
    {
        for (int j = 0; j <= jmax + 1; j++)
        {
            pic[i][j] = (((1 << CENTER) & Flag[i][j]) == 0) *
                        FLUID_PIXEL // if the cell is fluid then set pic value to FLUID_PIXEL
                        + (((1 << NSBIT) & Flag[i][j]) != 0) * NOSLIP_PIXEL // ...
                        + (((1 << FSBIT) & Flag[i][j]) != 0) * FREESLIP_PIXEL // ...
                        + (((1 << OFBIT) & Flag[i][j]) != 0) * OUTFLOW_PIXEL // ...
                        + (((1 << IFBIT) & Flag[i][j]) != 0) * INFLOW_PIXEL // ...
                        + (((1 << CBIT) & Flag[i][j]) != 0) *
                          COUPLING_PIXEL; // if the cell is coupling then set pic value to COUPLING_PIXEXL
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
    char *buffer = (char *) malloc((xsize * 2 + 1) * sizeof(char));
    
    for (int j = ysize - 1; j >= 0; j--)
    {
        int i = 0;
        
        for (i = 0; i <= 2 * xsize - 1; i = i + 2)
        {
            // fill in line buffer including spaces inbetween the numbers
            int value = pgm[i / 2][j];
            int character = 48 + value;
            buffer[i] = character;
            buffer[i + 1] = ' ';
        }
        
        // add newline character
        buffer[2 * xsize - 1] = '\n';
        // add end of string character
        buffer[2 * xsize] = '\0';
        
        fprintf(fp, "%s", buffer);
        //printf("%s\n",buffer);
    }
    
    // on Windows sth strange here - infinite loop (?) how to correctly free the buffer (?)
    free(buffer);
    buffer = NULL;
    fclose(fp);
    //printf("EOF");
}


int flipToFluid(double **U, double **V, int **Flag, int i, int j, int *obstacleBudget)
{
    if (isFluid(Flag[i][j]))
    {
        return 0;
    }
    Flag[i][j] = 0
                 + (1 << TOP) * isObstacle(Flag[i][j + 1])
                 + (1 << BOT) * isObstacle(Flag[i][j - 1])
                 + (1 << LEFT) * isObstacle(Flag[i - 1][j])
                 + (1 << RIGHT) * isObstacle(Flag[i + 1][j])
                 + (1 << TBIT) * (NEUMANN); // NEUMANN by default
    
    Flag[i + 1][j] -= (1 << LEFT);
    Flag[i - 1][j] -= (1 << RIGHT);
    Flag[i][j - 1] -= (1 << TOP);
    Flag[i][j + 1] -= (1 << BOT);
    
    logMsg(DEBUG, "Flipping to fluid: i=%d,j=%d", i, j);
    
    if (U!=NULL && V!=NULL)
    {
        U[i][j] = 0;
        V[i][j] = 0;
    }
    // Updating obstacle budget
    if (obstacleBudget!=NULL)
    {
        ++(*obstacleBudget);
    }
    return 1;
}

int flipToSolid(double **U, double **V, double **P, int **Flag, int i, int j, int *obstacleBudget)
{
    if (isObstacle(Flag[i][j]))
    {
        return 0;
    }
    Flag[i][j] = (1 << CENTER)
                 + (1 << NSBIT)
                 + (1 << TOP) * isObstacle(Flag[i][j + 1])
                 + (1 << BOT) * isObstacle(Flag[i][j - 1])
                 + (1 << LEFT) * isObstacle(Flag[i - 1][j])
                 + (1 << RIGHT) * isObstacle(Flag[i + 1][j])
                 + (1 << TBIT) * (NEUMANN); // NEUMANN by default
    
    
    // additional if() condition - for correctes of Flags[][] update in expandVortexSeeds()
    if (!isNeighbourObstacle(Flag[i+1][j], LEFT))
    {
        Flag[i + 1][j] += (1 << LEFT);
    }
    if (!isNeighbourObstacle(Flag[i - 1][j], RIGHT))
    {
        Flag[i - 1][j] += (1 << RIGHT);
    }
    if (!isNeighbourObstacle(Flag[i][j - 1], TOP))
    {
        Flag[i][j - 1] += (1 << TOP);
    }
    if (!isNeighbourObstacle(Flag[i][j + 1], BOT))
    {
        Flag[i][j + 1] += (1 << BOT);
    }
    
    
    logMsg(DEBUG, "Flipping to solid: i=%d,j=%d", i, j);
    
    if (U!=NULL && V!=NULL && P!=NULL)
    {
        U[i][j] = 0;
        U[i - 1][j] = 0;
        V[i][j] = 0;
        V[i][j - 1] = 0;
        P[i][j] = 0;
    
        // Now we need to make sure we don't leave any spurious tangential velocity on geometries that have become inner ones
        int cell = Flag[i][j];
        if (isNeighbourObstacle(cell, LEFT))
        {
            V[i - 1][j] = 0;
        }
        if (isNeighbourObstacle(cell, RIGHT))
        {
            V[i + 1][j] = 0;
        }
        if (isNeighbourObstacle(cell, BOT))
        {
            U[i][j - 1] = 0;
        }
        if (isNeighbourObstacle(cell, TOP))
        {
            U[i][j + 1] = 0;
        }
    }
    // Updating obstacle budget
    if (obstacleBudget!=NULL)
    {
        --(*obstacleBudget);
    }
    return 1;
}

/*void flipToSolidVortex(double **U, double **V, double **P, int **Flag, int i, int j)
{
    Flag[i][j] = (1 << CENTER)
                 + (1 << FSBIT)
                 + (1 << TOP)
                 + (1 << BOT) * isObstacle(Flag[i][j - 1])
                 + (1 << LEFT) * isObstacle(Flag[i - 1][j])
                 + (1 << RIGHT)
                 + (1 << TBIT) * (NEUMANN); // NEUMANN by default


    Flag[i+1][j] = (1 << CENTER)
                + (1 << FSBIT)
                + (1 << TOP)
                + (1 << BOT) * isObstacle(Flag[i + 1][j - 1])
                + (1 << LEFT)
                + (1 << RIGHT) * isObstacle(Flag[i + 2][j])
                + (1 << TBIT) * (NEUMANN); // NEUMANN by default

    Flag[i][j + 1] = (1 << CENTER)
                   + (1 << FSBIT)
                   + (1 << TOP) * isObstacle(Flag[i][j + 2])
                   + (1 << BOT)
                   + (1 << LEFT) * isObstacle(Flag[i - 1][j + 1])
                   + (1 << RIGHT)
                   + (1 << TBIT) * (NEUMANN); // NEUMANN by default

    Flag[i + 1][j + 1] = (1 << CENTER)
                     + (1 << FSBIT)
                     + (1 << TOP) * isObstacle(Flag[i +1][j + 2])
                     + (1 << BOT)
                     + (1 << LEFT)
                     + (1 << RIGHT) * isObstacle(Flag[i + 2][j + 1])
                     + (1 << TBIT) * (NEUMANN); // NEUMANN by default


    // plus update the surroundings of the whole big cell
    Flag[i-1][j] += (1 << RIGHT)*isObstacle(Flag[i][j]);
    Flag[i-1][j + 1] += (1 << RIGHT)*isObstacle(Flag[i][j + 1]);
    Flag[i][j -1] += (1 << TOP)*isObstacle(Flag[i][j]);
    Flag[i][j + 2] += (1 << BOT)*isObstacle(Flag[i][j + 1]);
    Flag[i+1][j - 1] += (1 << TOP)*isObstacle(Flag[i +1][j]);
    Flag[i+1][j + 2] += (1 << BOT)*isObstacle(Flag[i + 1][j +1]);
    Flag[i + 2][j] += (1 << LEFT)*isObstacle(Flag[i +1][j]);
    Flag[i + 2][j + 1] += (1 << LEFT)*isObstacle(Flag[i+1][j+1]);

    logMsg(DEBUG,"Flipping to solid: i=%d,j=%d",i,j);
    logMsg(DEBUG,"Flipping to solid: i=%d,j=%d",i+1,j);
    logMsg(DEBUG,"Flipping to solid: i=%d,j=%d",i,j+1);
    logMsg(DEBUG,"Flipping to solid: i=%d,j=%d",i+1,j+1);

    U[i][j] = 0;
    V[i][j] = 0;
    P[i][j] = 0;
    U[i+1][j] = 0;
    V[i+1][j] = 0;
    P[i+1][j] = 0;
    U[i][j+1] = 0;
    V[i][j+1] = 0;
    P[i][j+1] = 0;
    U[i+1][j+1] = 0;
    V[i+1][j+1] = 0;
    P[i+1][j+1] = 0;
}*/

int isVelocityMagnitudeBelowThreshold(double eps, double leftVelocity, double rightVelocity, double bottomVelocity,
                                      double topVelocity)
{
    if (sqrt(pow((leftVelocity + rightVelocity) * 0.5, 2) + pow((bottomVelocity + topVelocity) * 0.5, 2)) < eps)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

int checkPressure(double percent, double **P, int Flag, int i, int j)
{
    double yaxis = isNeighbourFluid(Flag, TOP) * P[i][j + 1] + isNeighbourFluid(Flag, BOT) * P[i][j - 1];
    double xaxis = isNeighbourFluid(Flag, RIGHT) * P[i + 1][j] + isNeighbourFluid(Flag, LEFT) * P[i - 1][j];
    
    if (fabs((yaxis - xaxis) / fmax(xaxis, yaxis)) > percent)
    {
        return 1;
    }
    
    return 0;
}

double getRelativeVelocityThreshold(double maxU, double maxV, double percent)
{
    return sqrt(pow(maxU, 2) + pow(maxV, 2)) * percent;
}

int isVelocityGreaterThanRelativeThreshold(double velocity, double maxU, double maxV, double percent)
{
    double velocityThreshold = getRelativeVelocityThreshold(maxU, maxV, percent);
    if (velocity > velocityThreshold)
    {
        return 1;
    }
    return 0;
}

int isVelocityAboveRelativeThreshold(int isFlip, double percent, double **U, double **V, int Flag, int i, int j,
                                     double maxU, double maxV)
{
    if (isFlip)
    {
        return 1;
    }
    double left_velocity = sqrt(pow(U[i - 1][j], 2) + pow(V[i - 1][j], 2)) * isNeighbourFluid(Flag, LEFT);
    double right_velocity = sqrt(pow(U[i + 1][j], 2) + pow(V[i + 1][j], 2)) * isNeighbourFluid(Flag, RIGHT);
    double top_velocity = sqrt(pow(U[i][j + 1], 2) + pow(V[i][j + 1], 2)) * isNeighbourFluid(Flag, TOP);
    double bottom_velocity = sqrt(pow(U[i][j - 1], 2) + pow(V[i][j - 1], 2)) * isNeighbourFluid(Flag, BOT);
    double velocity = left_velocity + right_velocity + top_velocity + bottom_velocity;
    
    if (isVelocityGreaterThanRelativeThreshold(velocity, maxU, maxV, percent))
    {
        return 1;
    }
    int val = 0;
    if (isCorner(Flag))
    {
        double yaxis_velocity = top_velocity + bottom_velocity;
        double xaxis_velocity = left_velocity + right_velocity;
        
        if (fabs((yaxis_velocity - xaxis_velocity) / fmin(yaxis_velocity, xaxis_velocity)) > percent)
        {
            val = 1;
        }
//        if (fabs(fmin(yaxis_velocity, xaxis_velocity) / fmax(yaxis_velocity, xaxis_velocity)) > 0.01) //test
//        {
//            val = 1;
//        }
    }
    return val;
}


bool isCellUpstreamOfObstacle(int i, int j, int **Flags, double **U, double **V, int imax, int jmax)
{
    if (i<2 || i>imax-1 || j<2 || j>jmax-1)
    {
        return 1;
    }
    int cell = Flags[i][j];
    if (!hasAtLeastOneObstacleNeighbour(cell))
    {
        return 0;
    }
    double u=0, v=0;
    if (isNeighbourObstacle(cell, RIGHT))
    {
        u = U[i-2][j];
        if (fabs(u) < fabs(V[i-2][j]))
            u = 0;
    }
    else if (isNeighbourObstacle(cell, LEFT))
    {
        u = -U[i+1][j]; // sign flip necessary for uniform conditions
        if (fabs(u) < fabs(V[i+1][j]))
            u = 0;
    }
    if (isNeighbourObstacle(cell, TOP))
    {
        v = V[i][j-2];
        if (fabs(v) < fabs(U[i][j-2]))
            v = 0;
    }
    else if (isNeighbourObstacle(cell, BOT))
    {
        v = -V[i][j+1]; // sign flip necessary for uniform conditions
        if (fabs(v) < fabs(U[i][j+1]))
            v = 0;
    }
    if (fabs(u) > fabs(v))
    {
        return u > 0;
    }
    else
    {
        return v > 0;
    }
}

void update_pgm(int imax, int jmax, int *noFluidCells, int **Flag, double **P, double **U, double **V, double minVelocity,
                double maxVelocity, double percentPressure, double percentVelocity, double maxU, double maxV, int k,
                int *obstacleBudget, double dx, double dy, int isPressure, int isVelocity, int isUpstreamCheckEnabled,
                double downstreamVelocityFactor)
{
    int addedCellsCounter = 0;
    int removedCellsCounter = 0;
    int cell;
    
    int istart, jstart, ibound, jbound, iflip, jflip, icounter, jcounter;
    
    icounter = k & 1;
    jcounter = (k >> 1) & 1;
    
    istart = (imax - 1) * icounter + 1;
    jstart = (jmax - 1) * jcounter + 1;
    
    ibound = (imax + 1) * !icounter;
    jbound = (jmax + 1) * !jcounter;
    
    iflip = pow(-1, icounter);
    jflip = pow(-1, jcounter);
    
    // Now allocate aux matrix to keep trace of just flipped geometries, to prevent cascade flipping!
    int **justFlipped = imatrix(0, imax + 1, 0, jmax + 1);

    
    // First pass on matrix is for switching to fluid, so we earn some budget
    for (int i = istart; i * iflip < ibound; i += iflip)
    {
        for (int j = jstart; j * jflip < jbound; j += jflip)
        {
            cell = Flag[i][j];
            if (isGeometryConstant(cell)) // If current cell cannot be changed, skip to next one
            {
                continue;
            }
            if (isObstacle(cell))
            {
                int isFlip = 0;
                if(isCorner(cell)){
                    if(isPressure){
                        isFlip = checkPressure(percentPressure, P, cell, i , j);
                    }
                    if(isVelocity){
                        isFlip += isVelocityAboveRelativeThreshold(isFlip, percentVelocity, U, V, cell, i , j, maxU, maxV);
                        isFlip += !isVelocityMagnitudeBelowThreshold(maxVelocity, U[i - 1][j], U[i][j], V[i][j - 1], V[i][j]); // EXPERIMENTAL!
                    }
                }
                if (isFlip)
                {
                    flipToFluid(U, V, Flag, i, j, obstacleBudget);
                    justFlipped[i][j] = 1;
                    (*noFluidCells)++;
                    ++removedCellsCounter;
                }
            }
            // This is done below
            else if ( // here we are fluid
                    (*obstacleBudget > 0)
                    && (
                            (doesNeighbourAllowFlipToSolid(Flag, i, j, RIGHT, imax, jmax) && !justFlipped[i + 1][j])
                            || (doesNeighbourAllowFlipToSolid(Flag, i, j, LEFT, imax, jmax) && !justFlipped[i - 1][j])
                            || (doesNeighbourAllowFlipToSolid(Flag, i, j, TOP, imax, jmax) && !justFlipped[i][j + 1])
                            || (doesNeighbourAllowFlipToSolid(Flag, i, j, BOT, imax, jmax) && !justFlipped[i][j - 1])
                    )
                    )
            {
                double thresholdVelocity = minVelocity;
                if(isUpstreamCheckEnabled && !isCellUpstreamOfObstacle(i, j, Flag, U, V, imax, jmax))
                {
                    thresholdVelocity = minVelocity/downstreamVelocityFactor;
                }
                int isFlip = isVelocityMagnitudeBelowThreshold(thresholdVelocity, U[i - 1][j], U[i][j], V[i][j - 1], V[i][j]);
//                double gradThreshold = 6;
//                isFlip += isGradientAtObstacleAboveThreshold(U, V, dx, dy, cell, i, j, gradThreshold);
                if (isFlip)
                {
                    flipToSolid(U, V, P, Flag, i, j, obstacleBudget);
                    justFlipped[i][j] = 1;
                    (*noFluidCells)--;
                    ++addedCellsCounter;
                }
            }
        }
    }
//    // Second pass on matrix is for switching to solid, so we can directly use the budget earned in previous pass
//    for (int i = istart; i * iflip < ibound; i += iflip)
//    {
//        for (int j = jstart; j * jflip < jbound; j += jflip)
//        {
//            cell = Flag[i][j];
//            if (isGeometryConstant(cell)) // If current cell cannot be changed, skip to next one
//            {
//                continue;
//            }
//            if (isFluid(cell)
//                    && (*obstacleBudget > 0)
//                    && ((doesNeighbourAllowFlipToSolid(Flag, i, j, RIGHT, imax, jmax) && !justFlipped[i + 1][j])
//                        || (doesNeighbourAllowFlipToSolid(Flag, i, j, LEFT, imax, jmax) && !justFlipped[i - 1][j])
//                        || (doesNeighbourAllowFlipToSolid(Flag, i, j, TOP, imax, jmax) && !justFlipped[i][j + 1])
//                        || (doesNeighbourAllowFlipToSolid(Flag, i, j, BOT, imax, jmax) && !justFlipped[i][j - 1])
//                    )
//                    )
//            {
//                int isFlip = isVelocityMagnitudeBelowThreshold(minVelocity, U[i - 1][j], U[i][j], V[i][j - 1], V[i][j]);
//                if (isFlip)
//                {
//                    flipToSolid(U, V, P, Flag, i, j, obstacleBudget);
//                    justFlipped[i][j] = 1;
//                    (*noFluidCells)--;
//                    ++addedCellsCounter;
//                }
//            }
//        }
//    }
    
    free_imatrix(justFlipped, 0, imax + 1, 0, jmax + 1);
    
    if (removedCellsCounter > 0 || addedCellsCounter > 0)
    {
        logMsg(PRODUCTION,"Updating PGM: flipped=%d, added=%d, removed=%d, balance=%d, obstacleBudget=%d, noFluidCells=%d",
               addedCellsCounter+removedCellsCounter,
               addedCellsCounter,
               removedCellsCounter,
               addedCellsCounter-removedCellsCounter,
               *obstacleBudget,
               *noFluidCells);
    }
}

int isGradientAtObstacleAboveThreshold(double **U, double **V, double dx, double dy, int cell, int i, int j,
                                       double gradThreshold)
{
    int val = 0;
    if (hasAtLeastOneObstacleNeighbour(cell))
    {
        double xGrad = fabs(U[i][j] - U[i-1][j]) / (dx*fmax(U[i][j], U[i-1][j]));
        double yGrad = fabs(V[i][j] - V[i][j-1]) / (dy*fmax(V[i][j], V[i][j-1]));
        double gradient = sqrt(xGrad * xGrad + yGrad * yGrad);
        if (gradient > gradThreshold)
        {
            logMsg(DEBUG, "Gradient at obstacle above threshold: grad=%f, threshold=%f", gradient, gradThreshold);
            val = 1;
        }
    }
    return val;
}

void expandVortexSeeds(int imax, int jmax, int *noFluidCells, double **U, double **V, double **P, int **Flag,
                       int vortexAreaThreshold, double swirlThresholdStrength, int *obstacleBudget)
{
    //Every point in the plane is checked in the same way
    //Except when a Vortex is encountered, i.e. Vortex != 0, then the furthermost interpolation point in this direction is taken
    
    int left_shift = 0;
    int right_shift = 0;
    int top_shift = 0;
    int bottom_shift = 0;
    
    double **Vortex = matrix(0, imax + 1, 0, jmax + 1);
    int **VortexAreaTags = imatrix(0, imax + 1, 0, jmax + 1);
    
    int noVortex = 0;
    int noVortexNew = 0;
    int noVortexOld = -1;
    
    //iterate until it "converges"
    while (noVortexNew != noVortexOld)
    {
        noVortex = 0;
        noVortexOld = noVortexNew;
        
        for (int i = 2; i < imax; i++)
        {
            for (int j = 2; j < jmax; j++)
            {
                left_shift = 1;
                right_shift = 1;
                top_shift = 1;
                bottom_shift = 1;
                
                int cell = Flag[i][j];
                if (isGeometryConstant(cell)) // If current cell cannot be changed, skip to next one
                {
                    continue;
                }
                
                if ((isObstacle(cell))
                    || (isNeighbourObstacle(cell, RIGHT))
                    || (isNeighbourObstacle(cell, TOP))
                        )
                {
                    //printf("Cell i=%d, j=%d ignored\n", i, j);
                    // ignore obstacle cells
                    // and fluid cells immediately neighbouring with obstacles on the right and top
                    // due to a "stencil" for vortex seed
                }
                else
                {
                    
                    if (Vortex[i][j] == 0)
                    {
                        //find rightmost seed
                        while (i + right_shift < imax-1
                               && (Vortex[i + right_shift][j] != 0)
                               && ( fsign(Vortex[i + right_shift][j]) == fsign(Vortex[i + right_shift + 1][j])))
                        {
                            right_shift++;
                            if (isObstacle(Flag[i + right_shift][j]))
                            {
                                break;
                            }
                        }
//                        right_shift--;
    
                        //find leftmost seed
                        while (i - left_shift > 2
                               && (Vortex[i - left_shift][j] != 0)
                               && (fsign(Vortex[i - left_shift][j]) == fsign(Vortex[i - left_shift - 1][j])))
                        {
                            left_shift++;
                            if (isObstacle(Flag[i - left_shift][j]))
                            {
                                break;
                            }
                        }
//                        left_shift--;
    
                        //find topmost seed
                        while (j + top_shift < jmax-1
                               && (Vortex[i][j + top_shift] != 0)
                               && (fsign(Vortex[i][j + top_shift]) == fsign(Vortex[i][j + top_shift + 1])))
                        {
                            top_shift++;
                            if (isObstacle(Flag[i][j+top_shift]))
                            {
                                break;
                            }
                        }
//                        top_shift--;
    
                        //find bottommost seed
                        while (j - bottom_shift > 2
                               && (Vortex[i][j - bottom_shift] != 0)
                               && (fsign(Vortex[i][j - bottom_shift]) == fsign(Vortex[i][j - bottom_shift - 1])))
                        {
                            bottom_shift++;
                            if (isObstacle(Flag[i][j-bottom_shift]))
                            {
                                break;
                            }
                        }
//                        bottom_shift--;
    
                            // compute approx vortex spanned area
                        int h_shift = left_shift + right_shift;
                        int v_shift = bottom_shift + top_shift;
                        int vortexArea = h_shift * v_shift;
//                        h_shift /= 2;
//                        v_shift /= 2;
//                        int min_shift = min(h_shift, v_shift);
                        int min_h_shift = min(left_shift, right_shift);
                        int min_v_shift = min(bottom_shift, top_shift);
                        int avg_h_shift = (left_shift + right_shift) / 2;
                        int avg_v_shift = (bottom_shift + top_shift) / 2;
                        int corrected_left_shift = min(left_shift, avg_h_shift);
                        int corrected_right_shift = min(right_shift, avg_h_shift);
                        int corrected_bottom_shift = min(bottom_shift, avg_v_shift);
                        int corrected_top_shift = min(top_shift, avg_v_shift);
    
                        // interpolate right point

//                        double u_right = (U[i + min_h_shift][j] + U[i + min_h_shift][j + 1]) * 0.5;
//                        double v_right = (V[i + min_h_shift][j] + V[(i + min_h_shift) + 1][j]) * 0.5;
//                        // interpolate left point
//                        double u_left = (U[i - min_h_shift][j] + U[i - min_h_shift][j + 1]) * 0.5;
//                        double v_left = (V[(i - min_h_shift) - 1][j] + V[(i - min_h_shift)][j]) * 0.5;
//                        // interpolate top point
//                        double u_top = (U[i][j + min_v_shift] + U[i][(j + min_v_shift) + 1]) * 0.5;
//                        double v_top = (V[i][j + min_v_shift] + V[i + 1][j + min_v_shift]) * 0.5;
//                        // interpolate bottom point
//                        double u_bottom = (U[i][(j - min_v_shift) - 1] + U[i][(j - min_v_shift)]) * 0.5;
//                        double v_bottom = (V[i][j - min_v_shift] + V[i + 1][j - min_v_shift]) * 0.5;
    
                        double u_right = (U[i + corrected_right_shift][j] + U[i + corrected_right_shift][j + 1]) * 0.5;
                        double v_right = (V[i + corrected_right_shift][j] + V[(i + corrected_right_shift) + 1][j]) * 0.5;
                        // interpolate left point
                        double u_left = (U[i - corrected_left_shift][j] + U[i - corrected_left_shift][j + 1]) * 0.5;
                        double v_left = (V[(i - corrected_left_shift) - 1][j] + V[(i - corrected_left_shift)][j]) * 0.5;
                        // interpolate top point
                        double u_top = (U[i][j + corrected_top_shift] + U[i][(j + corrected_top_shift) + 1]) * 0.5;
                        double v_top = (V[i][j + corrected_top_shift] + V[i + 1][j + corrected_top_shift]) * 0.5;
                        // interpolate bottom point
                        double u_bottom = (U[i][(j - corrected_bottom_shift) - 1] + U[i][(j - corrected_bottom_shift)]) * 0.5;
                        double v_bottom = (V[i][j - corrected_bottom_shift] + V[i + 1][j - corrected_bottom_shift]) * 0.5;
    
                        // calculate signum values
                        int sgn_v_left = 0;
                        int sgn_v_right = 0;
                        int sgn_u_top = 0;
                        int sgn_u_bottom = 0;
    
                        if (v_left > 0)
                        { sgn_v_left = 1; }
                        else if (v_left < 0) sgn_v_left = -1;
    
                        if (v_right > 0)
                        { sgn_v_right = 1; }
                        else if (v_right < 0) sgn_v_right = -1;
    
                        if (u_top > 0)
                        { sgn_u_top = 1; }
                        else if (u_top < 0) sgn_u_top = -1;
    
                        if (u_bottom > 0)
                        { sgn_u_bottom = 1; }
                        else if (u_bottom < 0) sgn_u_bottom = -1;
                        
                        
                        // check if it qualifies as a vortex seed
                        Vortex[i][j] = 0; // no rotation by default
                        if ((sgn_v_left + sgn_v_right + sgn_u_bottom + sgn_u_top == 0)
                            && (sgn_v_left + sgn_u_top != 0))
                        {
                            double strength = (fabs(v_left) + fabs(v_right) + fabs(u_top) + fabs(u_bottom))/4;
//                            if (strength > swirlThresholdStrength)
//                            {
                                noVortex++;
                                // determine in which direction it swirls
                                if (sgn_v_left == 1)
                                {
                                    Vortex[i][j] = -1 * strength; // right swirl
                                }
                                else
                                {
                                    Vortex[i][j] = 1 * strength; // left swirl
                                }
                                // Now tag the entire vortex with its area
                                double curVortex = Vortex[i][j];
                                int curVortexArea = vortexArea;
                                int iStart = max(i - avg_h_shift, 1);
                                int iEnd = min(i + avg_h_shift, imax);
                                int jStart = max(j - avg_v_shift, 1);
                                int jEnd = min(j + avg_v_shift, jmax);
                                // Get max vortex strength in modulus
                                for (int x = iStart; x <= iEnd; ++x)
                                {
                                    for (int y = jStart; y <= jEnd; ++y)
                                    {
                                        double curCellStrength = Vortex[x][y];
                                        int curCellArea = VortexAreaTags[x][y];
                                        if (fsign(curCellStrength) == fsign(curVortex))
                                        {
                                            if (fabs(curCellStrength) > fabs(curVortex))
                                            {
                                                curVortex = curCellStrength;
                                            }
                                            if (curCellArea > curVortexArea)
                                            {
                                                curVortexArea = curCellArea;
                                            }
                                        }
                                    }
                                }
                                // Copy max vortex strength and area into all cells
                                for (int x = iStart; x <= iEnd; ++x)
                                {
                                    for (int y = jStart; y <= jEnd; ++y)
                                    {
                                        if (fsign(Vortex[x][y]) == fsign(curVortex))
                                        {
                                            Vortex[x][y] = curVortex;
                                            VortexAreaTags[x][y] = curVortexArea;
                                        }
                                    }
                                }
//                            }
//                            else
//                            {
//                                logMsg(DEBUG, "Vortex too weak, ignoring it: strength=%f, swirlThresholdStrength=%f", strength, swirlThresholdStrength);
//                            }
                        }
                    }
                }
                
            }
            
        }
        
        noVortexNew = noVortexOld + noVortex;
    }
    
    // now fill in the geometries
    int filledCells = 0;
    for (int i = 1; i < imax; i++)
    {
        for (int j = 1; j < jmax; j++)
        {
            int curVortexArea = VortexAreaTags[i][j];
//            if (*obstacleBudget <= 0)
//            {
//                break;
//            }
            if (curVortexArea == 0)
            {
                continue;
            }
            if (curVortexArea < vortexAreaThreshold)
            {
                logMsg(DEBUG, "Found vortex cell but skipping as belonging to too small vortex: i=%d, j=%d, vortexArea=%d", i, j, curVortexArea);
                continue;
            }
//            if (curVortexArea > *obstacleBudget)
//            {
//                logMsg(DEBUG, "Found vortex cell but skipping as we don't have enough budget: i=%d, j=%d, vortexArea=%d, obstacleBudget=%d", i, j, curVortexArea, *obstacleBudget);
//                continue;
//            }
            double curVortexStrength = Vortex[i][j];
            if (fabs(curVortexStrength) < swirlThresholdStrength)
            {
                logMsg(DEBUG, "Vortex too weak, ignoring it: strength=%f, swirlThresholdStrength=%f", curVortexStrength, swirlThresholdStrength);
                continue;
            }
            if (Vortex[i][j] != 0
                && !isGeometryConstant(Flag[i][j])
                && !isGeometryConstant(Flag[i][j+1])
                && !isGeometryConstant(Flag[i+1][j])
                && !isGeometryConstant(Flag[i+1][j+1])
//                && hasAtLeastOneObstacleNeighbour(Flag[i][j]) // we want to flip only adjacent to other obstacles
                )
            {
                int noFlippedCells = 0;
                noFlippedCells += flipToSolid(U, V, P, Flag, i, j, obstacleBudget);
                noFlippedCells += flipToSolid(U, V, P, Flag, i + 1, j, obstacleBudget);
                noFlippedCells += flipToSolid(U, V, P, Flag, i, j + 1, obstacleBudget);
                noFlippedCells += flipToSolid(U, V, P, Flag, i + 1, j + 1, obstacleBudget);
                (*noFluidCells) -= noFlippedCells;
                filledCells += noFlippedCells;
                logMsg(DEBUG, "Filled vortex cell: i=%d, j=%d, vortexArea=%d, vortexStrength=%f", i, j, curVortexArea, curVortexStrength);
            }
            else
            {
                logMsg(DEBUG, "Vortex cell could not be flipped because in forbidden area: i=%d, j=%d, vortexArea=%d", i, j, curVortexArea);
            }
        }
    }
    
    free_matrix(Vortex, 0, imax + 1, 0, jmax + 1);
    free_imatrix(VortexAreaTags, 0, imax + 1, 0, jmax + 1);
    
    if (filledCells > 0)
    {
        logMsg(PRODUCTION, "Vortices found! Filled %d vortex cells.", filledCells);
    }
}

void geometryFix(double **U, double **V, double **P, int **Flag, int imax, int jmax, int *noFluidCells, int *obstacleBudget)
{
    int fixedCellsCounter = 0;
    for (int i = 1; i < imax + 1; i++)
    {
        for (int j = 1; j < jmax + 1; j++)
        {
            int cell = Flag[i][j];
            if (isObstacle(cell))
            {
                // Check that not both top and bottom are fluid
                if ((isFluid(Flag[i][j + 1]) && isFluid(Flag[i][j - 1])))
                {
                    // double top = sqrt(pow(U[i][j + 1], 2) + pow(V[i][j + 1], 2));
                    // double bot = sqrt(pow(U[i][j - 1], 2) + pow(V[i][j - 1], 2));
                    
                    // if (top > bot) {
                    //     flipToSolid(U, V, P, Flag, i, j - 1);
                    // } else {
                    //     flipToSolid(U, V, P, Flag, i, j + 1);
                    // }
                    // --(*noFluidCells);
                    flipToFluid(U, V, Flag, i, j, obstacleBudget);
                    ++(*noFluidCells);
                    ++fixedCellsCounter;
                }
                // Check that not both left and right are fluid
                if ((isFluid(Flag[i - 1][j]) && isFluid(Flag[i + 1][j])))
                {
                    // double left = sqrt(pow(U[i - 1][j], 2) + pow(V[i - 1][j], 2));
                    // double right = sqrt(pow(U[i + 1][j], 2) + pow(V[i + 1][j], 2));
                    // if (right > left) {
                    //     flipToSolid(U, V, P, Flag, i - 1, j);
                    // } else {
                    //     flipToSolid(U, V, P, Flag, i + 1, j);
                    // }
                    // --(*noFluidCells);
                    if (isObstacle(Flag[i][j]))
                    {
                        flipToFluid(U, V, Flag, i, j, obstacleBudget);
                        ++(*noFluidCells);
                        ++fixedCellsCounter;
                    }
                    
                }
            }
        }
    }
    
    for (int i = 1; i < imax + 1; i++)
    {
        for (int j = 1; j < jmax + 1; j++)
        {
            int cell = Flag[i][j];
            // This is to avoid 1-cell holes
            if (isFluid(cell)
                && isNeighbourObstacle(cell, LEFT)
                && isNeighbourObstacle(cell, RIGHT)
                && isNeighbourObstacle(cell, TOP)
                && isNeighbourObstacle(cell, BOT)
                    )
            {
                flipToSolid(U, V, P, Flag, i, j, obstacleBudget);
                --(*noFluidCells);
                ++fixedCellsCounter;
            }
        }
    }
    
    if (fixedCellsCounter > 0)
    {
        logMsg(PRODUCTION, "Forbidden geometries found! Fixed %d cells.", fixedCellsCounter);
    }
}

void velocityFix(double **U, double **V, int **Flag, int imax, int jmax)
{
    for (int i=1; i<=imax; ++i)
    {
        for (int j=1; j<=jmax; ++j)
        {
            if (isObstacle(Flag[i][j]))
            {
                U[i][j] = 0;
                U[i - 1][j] = 0;
                V[i][j] = 0;
                V[i][j - 1] = 0;
            }
        }
    }
}

void outputCalculation(double **U, double **V, int **Flags, int imax, int jmax, double *outflow)
{
    
    (*outflow) = 0;
    for (int i = 1; i <= imax; i++)
    {
        if (isOutflow(Flags[i][0]))
        {
            (*outflow) += sqrt(pow(U[i][0], 2) + pow(V[i][0], 2));
        }
        
        if (isOutflow(Flags[i][jmax + 1]))
        {
            (*outflow) += sqrt(pow(U[i][jmax + 1], 2) + pow(V[i][jmax + 1], 2));
        }
    }
    for (int j = 1; j <= jmax; j++)
    {
        if (isOutflow(Flags[0][j]))
        {
            (*outflow) += sqrt(pow(U[0][j], 2) + pow(V[0][j], 2));
        }
        
        if (isOutflow(Flags[imax + 1][j]))
        {
            (*outflow) += sqrt(pow(U[imax + 1][j], 2) + pow(V[imax + 1][j], 2));
        }
    }
}

void randomGeometryRemoval(int imax, int jmax, int *noFluidCells, int *obstacleBudget, int **Flags, double **P, double **U,
                           double **V, double removalProbability)
{
    for (int i=1; i<=imax; ++i)
    {
        for (int j=1; j<=jmax; ++j)
        {
            int cell = Flags[i][j];
            if (isObstacle(cell) && !isGeometryConstant(cell))
            {
                int r = rand();
                if (r < (int) round(RAND_MAX * removalProbability))
                {
                    int flipped = flipToFluid(U, V, Flags, i, j, obstacleBudget);
                    if (flipped)
                    {
                        (*noFluidCells)--;
                    }
                }
            }
        }
    }
}

