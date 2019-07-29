/* This file contains the functions to do the normal test. */
#include<math.h>
#include"functions.h"
#include"Datastruct.h"
#include<stdlib.h>
#include<stdio.h>

// input an array and its length and then calculate the JB value.
double JBTest(double *dist, int len)
{
	double mu;
	double K2, K3, K4;
	double Skew, Kurt, JB;
	int i;

	mu = 0.0;
	// calculate the mean
	for (i = 0; i < len; i++)
		mu += dist[i];

	mu /= len;

	K2 = 0.0;
	K3 = 0.0;
	K4 = 0.0;
	// Calculate K2, K3 and K4
	for (i = 0; i < len; i++)
	{
		K2 += (dist[i] - mu) * (dist[i] - mu);
		K3 += (dist[i] - mu) * (dist[i] - mu) * (dist[i] - mu);
		K4 += (dist[i] - mu) * (dist[i] - mu) * (dist[i] - mu) * (dist[i] - mu);
	}
	// Calculate Kurtosis and Skewness.
	Skew = K3 / pow(K2, 1.5);
	Kurt = K4 / (K2 * K2) - 3;

	// Calculate the JB value
	JB = Skew * Skew * 6.0 / len + Kurt * Kurt * 24.0 / len;

	return JB;
}

// Generate an assigned length random number array
// do the subtraction and scaling.
// conduct JBTest
int GenRandGArray(double *arr, int len, int Natom)
{
	int i, ki, yes;
	double mu, jb;
	double threshold;
	double v2, factor;

	// allocate memory
	if (arr==NULL)
		arr = (double*)malloc(len*sizeof(double));

	threshold = 0.1;
	yes = 1;

	while (yes)
	{
		mu = 0.0;
		// generate the  random number array
		for (i = 0; i < len; i++)
		{
			arr[i] = RandG(0,1);
			mu += arr[i];
		}
		mu /= len;

		v2 = 0.0;
		// subtraction
		for (i = 0; i < len; i++)
		{
			arr[i] -= mu;
			v2 += arr[i] * arr[i];
		}

		factor = sqrt(Natom / v2);
		// scaling
		for (i = 0; i < len; i++)
			arr[i] *= factor;

		jb = JBTest(arr,len);
		if (jb < threshold)
		{
			yes = 0;
			printf("Length: %d, Jarque-Bera Test Value: %.3e\n", len, jb);
		}
	}
	return 0;
}

// Write a function to continuously generate random number array until they pass the JBTest
int InitVelocityS(Sim *sim, Mat rho, Field *field, Coeff *coeff)
{
	int i, j, k, ci;
	int Nx, Ny, Nz;
	int lenx, leny, lenz;
	double *vs;
	FILE *fp;

	// retrive info
	Nx = rho.Nx;
	Ny = rho.Ny;
	Nz = rho.Nz;
	sim->Natom = 0;

	lenx = 0;
	leny = 0;
	lenz = 0;
	// first count the number for each direction.
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
			{
				if (rho.e[i][j][k] != 0)
					sim->Natom++;

				if (coeff->buox.e[i][j][k] == 1)
					lenx++;

				if (coeff->buoy.e[i][j][k] == 1)
					leny++;
				if (coeff->buoz.e[i][j][k] == 1)
					lenz++;
			}

	// allocate memory
	vs = (double*)malloc(lenx*sizeof(double));
	GenRandGArray(vs, lenx, sim->Natom);
	ci = 0;
	fp = fopen("Vx.log","a+");
	// Assign velocities to the simulation domain and record the initial values
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
			{
				if (coeff->buox.e[i][j][k] == 1)
				{
					fprintf(fp,"%e\n", vs[ci]);
					field->Vx.e[i][j][k] = vs[ci];
					ci++;
				}
			}
	free(vs);
	fclose(fp);

	vs = (double*)malloc(leny*sizeof(double));
	GenRandGArray(vs, leny, sim->Natom);
	ci = 0;
	fp = fopen("Vy.log","a+");
	// Assign velocities to the simulation domain
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
			{
				if (coeff->buoy.e[i][j][k] == 1)
				{
					fprintf(fp,"%e\n", vs[ci]);
					field->Vy.e[i][j][k] = vs[ci];
					ci++;
				}
			}
	free(vs);
	fclose(fp);

	vs = (double*)malloc(lenz*sizeof(double));
	GenRandGArray(vs, lenz, sim->Natom);
	ci = 0;
	fp = fopen("Vz.log","a+");
	// Assign velocities to the simulation domain
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
			{
				if (coeff->buoz.e[i][j][k] == 1)
				{
					fprintf(fp,"%e\n", vs[ci]);
					field->Vz.e[i][j][k] = vs[ci];
					ci++;
				}
			}
	free(vs);
	fclose(fp);

	return 0;
}
