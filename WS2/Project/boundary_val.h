#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__

/*
 * Auxiliary data structures to handle the boundary values.
 */
typedef enum BoundarySide
{
    T,
    B,
    L,
    R
} BoundarySide;

typedef enum BoundaryType
{
    DIRICHLET,
    INTERPOLATION,
    NEUMANN
} BoundaryType;

typedef struct BoundaryInfo
{
    BoundaryType typeU;
    BoundaryType typeV;
    char constU; // 1 means the value to apply is uniform.
    char constV; // 1 means the value to apply is uniform.
    double* valuesU;
    double* valuesV;
} BoundaryInfo;

/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(int imax, int jmax, double **U, double **V, int **Flags, BoundaryInfo boundaryInfo[4]);

#endif
