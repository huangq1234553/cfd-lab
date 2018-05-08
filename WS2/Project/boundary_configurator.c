//
// Created by tommaso on 07/05/18.
//

#include <memory.h>
#include "init.h"
#include "boundary_configurator.h"
#include "helper.h"

void configureBoundary(BoundaryInfo *boundaryInfo, BoundarySide boundarySide, const char *boundaryTypeStr,
                       const char *boundaryTempTypeStr, double u, double v, double temp, double qN, double k, double h)
{
    // Velocity boundary values
    if (strcmp(boundaryTypeStr, "NOSLIP") == 0)
    {
        initBoundaryInfo(boundaryInfo + boundarySide, DIRICHLET, DIRICHLET, 1, 1);
        *(boundaryInfo[boundarySide].valuesU) = 0;
        *(boundaryInfo[boundarySide].valuesV) = 0;
    }
    else if (strcmp(boundaryTypeStr, "MOVINGWALL") == 0)
    {
        initBoundaryInfo(boundaryInfo + boundarySide, DIRICHLET, DIRICHLET, 1, 1);
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
                initBoundaryInfo(boundaryInfo + boundarySide, DIRICHLET, NEUMANN, 1, 1);
                *(boundaryInfo[boundarySide].valuesU) = 0;
            case TOPBOUNDARY:
            case BOTTOMBOUNDARY:
                initBoundaryInfo(boundaryInfo + boundarySide, NEUMANN, DIRICHLET, 1, 1);
                *(boundaryInfo[boundarySide].valuesV) = 0;
        }
    }
    else if (strcmp(boundaryTypeStr, "INFLOW") == 0)
    {
        initBoundaryInfo(boundaryInfo + boundarySide, DIRICHLET, DIRICHLET, 1, 1);
        *(boundaryInfo[boundarySide].valuesU) = u;
        *(boundaryInfo[boundarySide].valuesV) = v;
    }
    else if (strcmp(boundaryTypeStr, "OUTFLOW") == 0)
    {
        initBoundaryInfo(boundaryInfo + boundarySide, NEUMANN, NEUMANN, 1, 1);
    }
    else
    {
        ERROR("Invalid velocity boundary type!");
    }
    
    // Temperature boundary values
    if (strcmp(boundaryTempTypeStr, "DIRICHLET") == 0)
    {
        boundaryInfo[boundarySide].typeT = DIRICHLET;
        *(boundaryInfo[boundarySide].valuesT) = temp;
    }
    else if (strcmp(boundaryTempTypeStr, "NEUMANN") == 0)
    {
        boundaryInfo[boundarySide].typeT = NEUMANN;
        boundaryInfo[boundarySide].coeff = h*qN/k;
    }
    else
    {
        ERROR("Invalid temperature boundary type!");
    }
}
