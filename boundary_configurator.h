//
// Created by tommaso on 07/05/18.
//

#ifndef SIM_BOUNDARY_CONFIGURATOR_H
#define SIM_BOUNDARY_CONFIGURATOR_H

typedef enum GeometryPixelValue
{
    NOSLIP_PIXEL = 0,
    FREESLIP_PIXEL = 1,
    OUTFLOW_PIXEL = 2,
    INFLOW_PIXEL = 3,
    COUPLING_PIXEL = 4,
    FLUID_PIXEL = 6,
} GeometryPixelValue;

void mapVelocityPixelValueToBoundaryTypeStr(GeometryPixelValue pixel, char* boundaryTypeStr);

void getVelocityBoundaryTypesFromExtendedGeometryFile(int** geometry, int imax, int jmax,
                                                      char* leftBoundaryType, char* rightBoundaryType,
                                                      char* topBoundaryType, char* bottomBoundaryType);

void configureVelocityBoundary(BoundaryInfo *boundaryInfo, const BoundarySide boundarySide,
                               const char *boundaryTypeStr,
                               const double u, const double v);

void configureTemperatureBoundary(BoundaryInfo *boundaryInfo, const BoundarySide boundarySide,
                                  const char *boundaryTempTypeStr,
                                  const double temperature, const double qN, const double k, const double h);

void configureBoundary(BoundaryInfo *boundaryInfo, const BoundarySide boundarySide, const char *boundaryTypeStr,
                       const char *boundaryTempTypeStr, double u, double v, double temp, double qN, double k,
                       double h);

#endif //SIM_BOUNDARY_CONFIGURATOR_H
