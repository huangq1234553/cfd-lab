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

// To be removed down the road
//typedef enum BoundaryType
//{
//    DIRICHLET,
//    NEUMANN
//} BoundaryType;
//
//typedef enum HLBoundaryType
//{
//    NOSLIP,
//    FREESLIP,
//    MOVINGWALL,
//    INFLOW,
//    OUTFLOW
//} HLBoundaryType;

typedef struct BoundaryInfo
{
    double *valuesDirichletU;
    double *valuesDirichletV;
    double *valuesDirichletT;
    double coeff; // This is dx*qN/k or dy*qn/k depending on the boundary
} BoundaryInfo;

//// Initialize a BoundaryInfo object
//void initBoundaryInfo(BoundaryInfo *boundaryInfo, HLBoundaryType hlBoundaryType, BoundaryType typeU, BoundaryType typeV,
//                      int numValuesU, int numValuesV);

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
