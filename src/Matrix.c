// This file helps with definition of the matrix data structure and operations.
#include"functions.h"
#include"Datastruct.h"
#include<stdio.h>
#include<stdlib.h>

// Set the size of the matrix
int SetMatSize(Mat *M, int Nx, int Ny, int Nz)
{
	M->Nx = Nx;
	M->Ny = Ny;
	M->Nz = Nz;
	// normal exit
	return 0;
}

// Allocate memory for matrix and initialized as zeros or ones.
int MallocMat(Mat *M, double InitVal)
{
	int i, j, k;
	int Nx, Ny, Nz;

	// retrive size
	Nx = M->Nx;
	Ny = M->Ny;
	Nz = M->Nz;

	M->e = (double***)malloc(Nx*sizeof(double**));
	for (i = 0; i < Nx; i++)
	{
		M->e[i] = (double**)malloc(Ny*sizeof(double*));
		for (j = 0; j < Ny; j++)
		{
			M->e[i][j] = (double*)malloc(Nz*sizeof(double));
			for (k = 0; k < Nz; k++)
				M->e[i][j][k] = InitVal;
		}
	}
	// normal exit
	return 0;
}
