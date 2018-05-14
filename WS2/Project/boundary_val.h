#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__

/*
 * Auxiliary data structures to handle the boundary values.
 */
typedef enum BoundarySide
{
    TOPBOUNDARY,
    BOTTOMBOUNDARY,
    LEFTBOUNDARY,
    RIGHTBOUNDARY
} BoundarySide;

typedef enum BoundaryType
{
    DIRICHLET,
    NEUMANN
} BoundaryType;

typedef struct BoundaryInfo
{
    BoundaryType typeU;
    BoundaryType typeV;
    BoundaryType typeT;
    char constU; // 1 means the value to apply is uniform.
    char constV; // 1 means the value to apply is uniform.
    char constT; // 1 means the value to apply is uniform.
    double *valuesU;
    double *valuesV;
    double *valuesT;
    double coeff; // This is dx*qN/k or dy*qn/k depending on the boundary
} BoundaryInfo;

// Initialize a BoundaryInfo object
void initBoundaryInfo(BoundaryInfo *boundaryInfo, BoundaryType typeU, BoundaryType typeV,
                      int numValuesU, int numValuesV);

void freeAllBoundaryInfo(BoundaryInfo boundaryInfo[4]);

/**
 * The boundary values of the problem are set.
 */

void boundaryval(int imax, int jmax, double **U, double **V, double **T, int **Flags,
                 BoundaryInfo *boundaryInfo);

void setLeftBoundaryValues(int imax, int jmax, double **U, double **V, double **T, int **Flags,
                           BoundaryInfo *boundaryInfo);

void setRightBoundaryValues(int imax, int jmax, double **U, double **V, double **T, int **Flags,
                            BoundaryInfo *boundaryInfo);

void setTopBoundaryValues(int imax, int jmax, double **U, double **V, double **T, int **Flags,
                          BoundaryInfo *boundaryInfo);

void setBottomBoundaryValues(int imax, int jmax, double **U, double **V, double **T, int **Flags,
                             BoundaryInfo *boundaryInfo);

void setEdgeBoundaryValues(int imax, int jmax, double **U, double **V, int **Flags, int i, int j);

void setCenterBoundaryValues(int imax, int jmax, double **Q, int **Flags, int i, int j);

#endif
