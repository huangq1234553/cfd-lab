#include "boundary_val.h"
#include "helper.h"

void boundaryvalues(int imax, int jmax, double **U, double **V, int **Flags)
{
	// TODO: Enhance this part to support configurable boundary conditions
	// Boundary values at the domain boundary
	for (int j = 1; j <= jmax; j++) {
		// left boundary - no-slip condition
		U[0][j] = 0;
		// right boundary - no-slip condition
		U[imax][j] = 0;
		// left boundary - no-slip condition - no v component on right wall of the cell - we need to enforce the condition v = 0 at the boundary by interpolation
		V[0][j] = -V[1][j];
		// right boundary - no-slip condition - no v component on right wall of the cell - we need to enforce the condition v = 0 at the boundary by interpolation
		V[imax+1][j] = -V[imax][j];
	}

	for (int i = 1; i <= imax; i++) {
		// lower boundary - no-slip condition
		V[i][0] = 0;
		// upper boundary - moving wall condition - set v = 0
		V[i][jmax] = 0;
		// lower boundary - no-slip condition - no u component on the horizontal wall of the cell - we need to enforce the condition u = 0 at the boundary by interpolation
		U[i][0] = -U[i][1];
		// upper boundary - moving wall condition - no u component on the horizontal wall of the cell - we need to enforce the condition u = 1 at the boundary by interpolation
		U[i][jmax+1] = 2.0 - U[i][jmax];
	}
    
    // Boundary values at geometries in the internal part of the domain
    for (int i = 1; i <= imax; ++i)
    {
        for (int j =1; j <= jmax; ++j)
        {
            int cell = Flags[i][j];
            if (isObstacle(cell))
            {
                // Compute v
                if (!skipV(cell))
                {
                    if (isNeighbourFluid(cell,TOP))
                    {
                        V[i][j] = 0;
                    }
                    else
                    {
                        V[i][j] = -V[i+isNeighbourObstacle(cell,LEFT)-isNeighbourObstacle(cell,RIGHT)][j];
                    }
                }
                // Compute u
                if (!skipU(cell))
                {
                    if (isNeighbourFluid(cell,RIGHT))
                    {
                        U[i][j] = 0;
                    }
                    else
                    {
                        U[i][j] = -U[i][j+isNeighbourObstacle(cell,BOT)-isNeighbourObstacle(cell,TOP)];
                    }
                }
            }
            else // if (isFluid(cell))
            {
                //compute V
                if (isNeighbourObstacle(cell,TOP))
                {
                    V[i][j] = 0;
                }
                //compute U
                if (isNeighbourObstacle(cell,RIGHT))
                {
                    U[i][j] = 0;
                }
            }
        }
    }
}

//eof
