//
// Created by tommaso on 07/05/18.
//

#include <memory.h>
#include "init.h"
#include "boundary_configurator.h"
#include "helper.h"

void
configureVelocityBoundary(BoundaryInfo *boundaryInfo, const BoundarySide boundarySide, const char *boundaryTypeStr,
                          const double u, const double v)
{
    // Velocity boundary values
    if (strcmp(boundaryTypeStr, "NOSLIP") == 0)
    {
        initBoundaryInfo(boundaryInfo + boundarySide, NOSLIP, DIRICHLET, DIRICHLET, 1, 1);
        *(boundaryInfo[boundarySide].valuesU) = 0;
        *(boundaryInfo[boundarySide].valuesV) = 0;
    }
    else if (strcmp(boundaryTypeStr, "MOVINGWALL") == 0)
    {
        initBoundaryInfo(boundaryInfo + boundarySide, MOVINGWALL, DIRICHLET, DIRICHLET, 1, 1);
        switch (boundarySide)
        {
            case LEFTBOUNDARY:
            case RIGHTBOUNDARY:
                *(boundaryInfo[boundarySide].valuesU) = 0;
                *(boundaryInfo[boundarySide].valuesV) = v;
            case TOPBOUNDARY:
            case BOTTOMBOUNDARY:
                *(boundaryInfo[boundarySide].valuesU) = u;
                *(boundaryInfo[boundarySide].valuesV) = 0;
        }
    }
    else if (strcmp(boundaryTypeStr, "FREESLIP") == 0)
    {
        switch (boundarySide)
        {
            case LEFTBOUNDARY:
            case RIGHTBOUNDARY:
                initBoundaryInfo(boundaryInfo + boundarySide, FREESLIP, DIRICHLET, NEUMANN, 1, 1);
                *(boundaryInfo[boundarySide].valuesU) = 0;
            case TOPBOUNDARY:
            case BOTTOMBOUNDARY:
                initBoundaryInfo(boundaryInfo + boundarySide, FREESLIP, NEUMANN, DIRICHLET, 1, 1);
                *(boundaryInfo[boundarySide].valuesV) = 0;
        }
    }
    else if (strcmp(boundaryTypeStr, "INFLOW") == 0)
    {
        initBoundaryInfo(boundaryInfo + boundarySide, INFLOW, DIRICHLET, DIRICHLET, 1, 1);
        *(boundaryInfo[boundarySide].valuesU) = u;
        *(boundaryInfo[boundarySide].valuesV) = v;
    }
    else if (strcmp(boundaryTypeStr, "OUTFLOW") == 0)
    {
        initBoundaryInfo(boundaryInfo + boundarySide, OUTFLOW, NEUMANN, NEUMANN, 1, 1);
    }
    else
    {
        THROW_ERROR("Invalid velocity boundary type!");
    }
}

void getVelocityBoundaryTypesFromExtendedGeometryFile(int **geometry, int imax, int jmax, char *leftBoundaryType,
                                                      char *rightBoundaryType, char *topBoundaryType,
                                                      char *bottomBoundaryType)
{
    // Here we need to loop on the entire boundaries because we never know if there is one spurious invalid one [tbc]
    for (int i = 1; i < imax+1; ++i)
    {
        mapVelocityPixelValueToBoundaryTypeStr((GeometryPixelValue) geometry[i][0], bottomBoundaryType);
        mapVelocityPixelValueToBoundaryTypeStr((GeometryPixelValue) geometry[i][jmax+1], topBoundaryType);
    }
    for (int j = 1; j < jmax+1; ++j)
    {
        mapVelocityPixelValueToBoundaryTypeStr((GeometryPixelValue) geometry[0][j], leftBoundaryType);
        mapVelocityPixelValueToBoundaryTypeStr((GeometryPixelValue) geometry[imax+1][j], rightBoundaryType);
    }
}

void mapVelocityPixelValueToBoundaryTypeStr(GeometryPixelValue pixel, char* boundaryTypeStr)
{
    switch (pixel)
    {
        case NOSLIP_PIXEL:
            strcpy(boundaryTypeStr, "NOSLIP");
            break;
        case FREESLIP_PIXEL:
            strcpy(boundaryTypeStr, "FREESLIP");
            break;
        case INFLOW_PIXEL:
            strcpy(boundaryTypeStr, "INFLOW");
            break;
        case OUTFLOW_PIXEL:
            strcpy(boundaryTypeStr, "OUTFLOW");
            break;
        default:
        {} // Else skip
    }
}

void configureTemperatureBoundary(BoundaryInfo *boundaryInfo, const BoundarySide boundarySide,
                                  const char *boundaryTempTypeStr, const double temperature, const double qN,
                                  const double k, const double h)
{
    // Temperature boundary values
    if (strcmp(boundaryTempTypeStr, "DIRICHLET") == 0)
    {
        boundaryInfo[boundarySide].typeT = DIRICHLET;
        *(boundaryInfo[boundarySide].valuesT) = temperature;
    }
    else if (strcmp(boundaryTempTypeStr, "NEUMANN") == 0)
    {
        boundaryInfo[boundarySide].typeT = NEUMANN;
        boundaryInfo[boundarySide].coeff = h * qN / k;
    }
    else
    {
        THROW_ERROR("Invalid temperature boundary type!");
    }
}

void configureBoundary(BoundaryInfo *boundaryInfo, const BoundarySide boundarySide, const char *boundaryTypeStr,
                       const char *boundaryTempTypeStr, double u, double v, double temp, double qN, double k,
                       double h)
{
    configureVelocityBoundary(boundaryInfo, boundarySide, boundaryTypeStr, u, v);
    configureTemperatureBoundary(boundaryInfo, boundarySide, boundaryTempTypeStr, temp, qN, k, h);
    
}
