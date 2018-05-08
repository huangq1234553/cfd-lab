//
// Created by tommaso on 07/05/18.
//

#ifndef SIM_BOUNDARY_CONFIGURATOR_H
#define SIM_BOUNDARY_CONFIGURATOR_H

void configureBoundary(BoundaryInfo *boundaryInfo, BoundarySide boundarySide, const char *boundaryTypeStr,
                       const char *boundaryTempTypeStr, double u, double v, double temp, double qN, double k, double h);

#endif //SIM_BOUNDARY_CONFIGURATOR_H
