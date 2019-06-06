// This file reads the density matrix and conduct parameter averaging algorithm for the boundary conditions, STILL NEEDS DOUBLE CHECK!!

#include<stdio.h>
#include<stdlib.h>
#include"Datastruct.h"
#include"functions.h"

// read in the density matrix
int ReadDensityMat(Mat *rho)
{
	FILE *fp;
	int i, j, k;

	fp = fopen("density.in","a+");
	// scanf the size of the matrix.
	fscanf(fp,"%d", &rho->Nx);
	fscanf(fp,"%d", &rho->Ny);
	fscanf(fp,"%d", &rho->Nz);

	// allocate memory for the density matrix and initialized with zeros
	MallocMat(rho, 0.0);

	// read in the density matrix
	for (k = 0; k < rho->Nz; k++)
		for (i = 0; i < rho->Nx; i++)
			for (j = 0; j < rho->Ny; j++)
				fscanf(fp,"%lf", &rho->e[i][j][k]);

	fclose(fp);
	// output the geometry
	fp = fopen("DenstiyMat.log","a+");
	for (k = 0; k < rho->Nz; k++)
	{
		for (i = 0; i < rho->Nx; i++)
		{
			for (j = 0; j < rho->Ny; j++)
				fprintf(fp,"%e  ", rho->e[i][j][k]);
			fprintf(fp,"\n");
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
	// normal exit
	return 0;
}
