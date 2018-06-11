#include <stdlib.h>
#include "precice_adapter.h"
#include "adapters/c/SolverInterfaceC.h"
#include "helper.h"
#include "boundary_val.h"

void write_checkpoint(double t, double **U, double **V, double **T, double **U_cp, double **V_cp,
                      double **T_cp, int imax, int jmax){
    for(int i=0; i < imax+1; i++){
        for(int j=0; j < jmax+1; j++){
            T_cp[i][j] = T[i][j];
            U_cp[i][j] = U[i][j];
            V_cp[i][j] = V[i][j];
        }
    }

}

void restore_checkpoint(double t, double **U, double **V, double **T, double **U_cp,
                        double **V_cp,
                        double **T_cp, int imax, int jmax){
    for(int i=0; i < imax+1; i++){
        for(int j=0; j < jmax+1; j++){
            T[i][j] = T_cp[i][j];
            U[i][j] = U_cp[i][j];
            V[i][j] = V_cp[i][j];
        }
    }
}


int *precice_set_interface_vertices(int imax, int jmax, double dx, double dy, double x_origin, double y_origin,
                                    int num_coupling_cells, int meshID, int **Flags, int dim){
    double* vertices = (double*) malloc(sizeof(double) * num_coupling_cells * dim);
    int* vertexIDs = (int*) malloc(sizeof(int) * num_coupling_cells);
    int count = 0;
    // Iterate over left boundary
    for(int j = 1; j <= jmax; j++){
        if(isCoupling(Flags[0][j])){
            vertices[count] =   x_origin + (-0.5)*dx;
            vertices[count+1] = y_origin + (j-0.5)*dy;
            vertices[count+2] = 0;
            count += 3;
        }
    }
    // Iterate over right boundary
    for(int j = 1; j <= jmax; j++){
        if(isCoupling(Flags[imax+1][j])){
            vertices[count] =   x_origin + (imax+1-0.5)*dx;
            vertices[count+1] = y_origin + (j-0.5)*dy;
            vertices[count+2] = 0;
            count += 3;
        }
    }
    // Iterate over top boundary
    for(int i = 1; i <= imax; i++){
        if(isCoupling(Flags[i][jmax+1])){
            vertices[count] =   x_origin + (i-0.5)*dx;
            vertices[count+1] = y_origin + (jmax+1-0.5)*dy;
            vertices[count+2] = 0;
            count += 3;
        }
    }
    // Iterate over bottom boundary
    for(int i = 1; i <= imax; i++){
        if(isCoupling(Flags[i][0])){
            vertices[count] =   x_origin + (i-0.5)*dx;
            vertices[count+1] = y_origin + (-0.5)*dy;
            vertices[count+2] = 0;
            count += 3;
        }
    }

    for(int i=1; i<=imax; i++){
        for(int j=1; j<=jmax; j++){
            if(isCoupling(Flags[i][j])){
                vertices[count] =   x_origin + (i-0.5)*dx;
                vertices[count+1] = y_origin + (j-0.5)*dy;
                vertices[count+2] = 0;
                count += 3;
            }
        }
    }

    precicec_setMeshVertices(meshID, num_coupling_cells, vertices, vertexIDs);
    free(vertices);

    return vertexIDs;

}

void precice_write_temperature(int imax, int jmax, int num_coupling_cells, double *temperature, int *vertexIDs,
                               int temperatureID, double **T, int **Flags){
    int count = 0;
    // Iterate over left boundary
    for(int j = 1; j <= jmax; j++){
        if(isCoupling(Flags[0][j])){
            temperature[count] = T[1][j];
            count += 1;
        }
    }
    // Iterate over right boundary
    for(int j = 1; j <= jmax; j++){
        if(isCoupling(Flags[imax+1][j])){
            temperature[count] = T[imax][j];
            count += 1;
        }
    }
    // Iterate over top boundary
    for(int i = 1; i <= imax; i++){
        if(isCoupling(Flags[i][jmax+1])){
            temperature[count] = T[i][jmax];
            count += 1;
        }
    }
    // Iterate over bottom boundary
    for(int i = 1; i <= imax; i++){
        if(isCoupling(Flags[i][0])){
            temperature[count] = T[i][1];
            count += 1;
        }
    }

    for(int i=1; i<=imax; i++){
        for(int j=1; j<=jmax; j++){
            if(isCoupling(Flags[i][j])){
                temperature[count] = T[i][j+1] * isNeighbourFluid(Flags[i][j], TOP)
                                        + T[i][j-1] * isNeighbourFluid(Flags[i][j], BOT);
                count += 1;
            }
        }
    }
    precicec_writeBlockScalarData(temperatureID, num_coupling_cells, vertexIDs, temperature);

}

void set_coupling_boundary(int imax, int jmax, double dx, double dy, double *heatflux, double **T, int **Flags){
    int count = 0;
    // Iterate over left boundary
    for(int j = 1; j <= jmax; j++){
        if(isCoupling(Flags[0][j])){
            T[0][j] = T[1][j] + dx*heatflux[count];
            count += 1;
        }
    }
    // Iterate over right boundary
    for(int j = 1; j <= jmax; j++){
        if(isCoupling(Flags[imax+1][j])){
            T[imax+1][j] = T[imax][j] + dx*heatflux[count];
            count += 1;
        }
    }
    // Iterate over top boundary
    for(int i = 1; i <= imax; i++){
        if(isCoupling(Flags[i][jmax+1])){
            T[i][jmax+1] = T[i][jmax] + dy*heatflux[count];
            count += 1;
        }
    }
    // Iterate over bottom boundary
    for(int i = 1; i <= imax; i++){
        if(isCoupling(Flags[i][0])){
            T[i][0] = T[i][1] + dy*heatflux[count];
            count += 1;
        }
    }

    for(int i=1; i<=imax; i++){
        for(int j=1; j<=jmax; j++){
            if(isCoupling(Flags[i][j])){
                T[i][j] = T[i][j+1] * isNeighbourFluid(Flags[i][j], TOP)
                            + T[i][j-1] * isNeighbourFluid(Flags[i][j], BOT)
                            + dy*heatflux[count];
                count += 1;
            }
        }
    }
}
