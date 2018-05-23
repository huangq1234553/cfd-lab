#include "boundary_val.h"

void boundaryvalues(int omg_i, int omg_j, int imax_local, int jmax_local, double **U, double **V) {
	
	if(omg_i == 0)
	{
		// left ghost layer
		// NOTE: U and V have different sizes, hence indexing shift
		for (int j = 1; j < jmax_local + 2; j++) {
			// left ghost layer - no-slip condition
			U[1][j] = 0;
			// left ghost layer - no-slip condition - no v component on right wall of the cell - we need to enforce the condition v = 0 at the boundary by interpolation
			V[0][j + 1] = -V[1][j + 1];
		}
	}
	else // if omg_i == iproc - 1
	{
		// right ghost layer
		// NOTE: U and V have different sizes, hence indexing shift
		for (int j = 1; j < jmax_local + 2; j++) {
			// right ghost layer - no-slip condition
			U[imax_local + 3 - 1][j] = 0;
			// right ghost layer - no-slip condition - no v component on right wall of the cell - we need to enforce the condition v = 0 at the boundary by interpolation
			V[imax_local + 2][j + 1] = -V[imax_local + 2 - 1][j + 1];
		}
	}

	if (omg_j == 0)
	{
		// bottom ghost layer
		// NOTE: U and V have different sizes, hence indexing shift
		for (int i = 1; i < imax_local + 2; i++) {
			// lower boundary - no-slip condition
			V[i][1] = 0;
			// lower boundary - no-slip condition - no u component on the horizontal wall of the cell - we need to enforce the condition u = 0 at the boundary by interpolation
			U[i + 1][0] = -U[i + 1][1];
		}
	}
	else // if omg_j = jproc - 1
	{
		// top ghost layer
		// NOTE: U and V have different sizes, hence indexing shift
		for (int i = 1; i < imax_local + 2; i++) {
			// upper boundary - moving wall condition - set v = 0
			V[i][jmax_local + 3 - 1] = 0;
			// upper boundary - moving wall condition - no u component on the horizontal wall of the cell - we need to enforce the condition u = 1 at the boundary by interpolation
			U[i + 1][jmax_local + 2] = 2.0 - U[i + 1][jmax_local + 2 - 1];
		}
	}
}


void boundaryvalues_FG(int omg_i, int omg_j, int imax_local, int jmax_local, double **F, double **G, double **U, double **V) {
	
	if(omg_i == 0)
	{
		// left ghost layer
		// NOTE: F and G have different sizes, hence indexing shift
		for (int j = 1; j < jmax_local + 2; j++) {
			// F left boundary condition
			F[1][j] = U[1][j];
		}
	}
	else // if omg_i == iproc - 1
	{
		// right ghost layer
		// NOTE: F and G have different sizes, hence indexing shift
		for (int j = 1; j < jmax_local + 2; j++) {
			// F right boundary condition
					//TODO: should the index not be [imax_local + 2][j]?
			F[imax_local+3][j] = U[imax_local+3][j];
		}
	}

	if (omg_j == 0)
	{
		// bottom ghost layer
		// NOTE: F and G have different sizes, hence indexing shift
		for (int i = 1; i < imax_local + 2; i++) {
			// G lower boundary condition
			G[i][1] = V[i][1];
		}
	}
	else // if omg_j = jproc - 1
	{
		// top ghost layer
		// NOTE: F and G have different sizes, hence indexing shift
		for (int i = 1; i < imax_local + 2; i++) {
			// G top boundary condition
			//TODO: should the index not be [i][jmax_local+2]?
			G[i][jmax_local+3] = V[i][jmax_local+3];
		}
	}
}

void boundaryvalues_P(int omg_i, int omg_j, int imax_local, int jmax_local, double **P)
{
    /* set boundary values */
    int i,j;

    if (omg_i == 0)
    {
        for(i = 1; i <= imax_local; i++) {
            P[i][0] = P[i][1];
        }
    }
    else
    {
        for(i = 1; i <= imax_local; i++) {
            P[i][jmax_local+1] = P[i][jmax_local];
        }
    }

    if(omg_j == 0){
        for(j = 1; j <= jmax_local; j++) {
            P[0][j] = P[1][j];
        }
    }
    else{
        for(j = 1; j <= jmax_local; j++) {
            P[imax_local+1][j] = P[imax_local][j];
        }
    }
}

//eof
