#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__

#include <stdbool.h>

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
    DIRICHLET=0,
    NEUMANN=1
} BoundaryType;

// To be removed down the road
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
    double *valuesDirichletP;
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

void setCenterBoundaryValues(int imax, int jmax, double **Q, int **Flags, int i, int j, bool isCoupled);

void setPressureOuterBoundaryValues(int imax, int jmax, double **P, int **Flags, const BoundaryInfo *boundaryInfo);

#endif
