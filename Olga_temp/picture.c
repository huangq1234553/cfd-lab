// picture.c : Defines the entry point for the console application.
//

#include <stdio.h>      /* printf, scanf, NULL */
#include <stdlib.h>     /* malloc, free, rand */

typedef enum Direction { CENTER = 0, TOP = 4, BOT = 3, LEFT = 2, RIGHT = 1 } Direction;
typedef enum BoundaryTypeBit { NSBIT = 5, FSBIT = 6, OFBIT = 7, IFBIT = 8, CBIT = 9, TBIT = 10 } BoundaryTypeBit;
typedef enum GeometryPixelValue { NOSLIP_PIXEL = 0, FREESLIP_PIXEL = 1, OUTFLOW_PIXEL = 2, INFLOW_PIXEL = 3, COUPLING_PIXEL = 4, FLUID_PIXEL = 6 } GeometryPixelValue;
typedef enum BoundaryType { DIRICHLET = 0, NEUMANN = 1 } BoundaryType;

void encode_flags(int imax, int jmax, int** Flag, int** pic)
{
	for (int i = 0; i <= imax; i++)
	{
		for (int j = 0; j <= jmax; j++)
		{
			Flag[i][j] = (1 << CENTER) * (pic[i][j] != FLUID_PIXEL)
				+ (1 << NSBIT) * (pic[i][j] == NOSLIP_PIXEL)
				+ (1 << FSBIT) * (pic[i][j] == FREESLIP_PIXEL)
				+ (1 << OFBIT) * (pic[i][j] == OUTFLOW_PIXEL)
				+ (1 << IFBIT) * (pic[i][j] == INFLOW_PIXEL)
				+ (1 << CBIT) * (pic[i][j] == COUPLING_PIXEL)
				+ (1 << TBIT) * (NEUMANN); // NEUMANN by default
		}
	}
}

void decode_flags(int imax, int jmax, int** Flag, int** pic)
{
	// decode flags and return gray colour map - go through the whole domain including boundaries, i.e. EXTENDED mode
	for (int i = 0; i <= imax; i++)
	{
		for (int j = 0; j <= jmax; j++)
		{
			pic[i][j] = (((1 << CENTER) & Flag[i][j]) == 0) * FLUID_PIXEL // if the cell is fluid then set pic value to FLUID_PIXEL
				+ (((1 << NSBIT) & Flag[i][j]) != 0) * NOSLIP_PIXEL // ...
				+ (((1 << FSBIT) & Flag[i][j]) != 0) * FREESLIP_PIXEL // ...
				+ (((1 << OFBIT) & Flag[i][j]) != 0) * OUTFLOW_PIXEL // ...
				+ (((1 << IFBIT) & Flag[i][j]) != 0) * INFLOW_PIXEL // ...
				+ (((1 << CBIT) & Flag[i][j]) != 0) * COUPLING_PIXEL; // if the cell is coupling then set pic value to COUPLING_PIXEXL
		}
	}
}

void setup_test_pgm(int** pgm, int imax, int jmax)
{
	for (int i = 0; i <= imax; i++)
	{
		for (int j = 0; j <= jmax; j++)
		{
			if (j == 0 && (i != 0 && i != imax))
				pgm[i][j] = INFLOW_PIXEL;
			else if (j == jmax && (i != 0 && i != imax))
				pgm[i][j] = OUTFLOW_PIXEL;
			else if (i == 0 || i == imax)
				pgm[i][j] = FREESLIP_PIXEL;
			else
				pgm[i][j] = FLUID_PIXEL;
		}
	}
}

void write_pgm(int xsize, int ysize, int** pgm, const char* filename)
{
	FILE *output = NULL;

	if ((output = fopen(filename, "wb")) != 0) {
		printf("write_pgm(): filename = %s created\n", filename);
	}
	else
	{
		printf("write_pgm(): unable to create filename = %s\n", filename);
	}

	fprintf(output, "P2\n");
	//printf("P2\n");

	fprintf(output, "# 0 - no-slip\n");
	//printf("# 0 - no-slip\n");

	fprintf(output, "# 1 - free-slip\n");
	//printf("# 1 - free-slip\n");

	fprintf(output, "# 2 - outflow\n");
	//printf("# 2 - outflow\n");

	fprintf(output, "# 3 - inflow\n");
	//printf("# 3 - inflow\n");

	fprintf(output, "# 4 - coupling\n");
	//printf("# 4 - coupling\n");

	fprintf(output, "# 6 - fluid\n");
	//printf("# 6 - fluid\n");

	fprintf(output, "%d %d\n", ysize, xsize);
	//printf("%d %d\n", xsize, ysize);

	fprintf(output, "%d\n", FLUID_PIXEL);
	//printf("%d\n", FLUID_PIXEL);

	// create buffer
	char* buffer = (char*)malloc((ysize * 2 + 1) * sizeof(char));

	for (int i = xsize - 1; i >= 0; i--)
	{
		int j = 0;

		for(j = 0; j <= 2 * ysize - 1; j = j+2)
		{
			// fill in line buffer including spaces inbetween
			buffer[j] = 48 + pgm[i][j/2];
			buffer[j + 1] = ' ';
		}

		buffer[2*ysize - 1] = '\n';
		buffer[2*ysize] = '\0';

		fprintf(output, "%s", buffer);
		//printf("%s\n",buffer);
	}

	free(buffer);
	buffer = NULL;

	fclose(output);
	// printf("EOF");
}

int main()
{
	// fix size of "pgm" inner domain
	int imax = 49;
	int jmax = 69;

	// allocate memory for matrices
	int** pgmM = (int**)malloc((imax + 1) * sizeof(int));
	int** flagM = (int**)malloc((imax + 1) * sizeof(int));
	int** pgm_retrievedM = (int**)malloc((imax + 1) * sizeof(int));

	for (int i = 0; i <= imax; i++)
	{
		pgmM[i] = (int*)malloc((jmax + 1) * sizeof(int));
		flagM[i] = (int*)malloc((jmax + 1) * sizeof(int));
		pgm_retrievedM[i] = (int*)malloc((jmax + 1) * sizeof(int));
	}

	// setup test "pgm"
	setup_test_pgm(pgmM, imax, jmax);

	// encode flags
	encode_flags(imax, jmax, flagM, pgmM);

	// decode flags
	decode_flags(imax, jmax, flagM, pgm_retrievedM);

	write_pgm(imax + 1, jmax + 1, pgmM, "test.pgm");
	write_pgm(imax + 1, jmax + 1, pgm_retrievedM, "retrieved.pgm");


	// free columns
	for (int i = 0; i <= imax; i++)
	{
		free(pgmM[i]);
		pgmM[i] = NULL;
		free(flagM[i]);
		flagM[i] = NULL,
		free(pgm_retrievedM[i]);
		pgm_retrievedM[i] = NULL;
	}

	// free rows
	free(pgmM);
	pgmM = NULL;
	free(pgm_retrievedM);
	pgm_retrievedM = NULL;
	free(flagM);
	flagM = NULL;

	return 0;
}

