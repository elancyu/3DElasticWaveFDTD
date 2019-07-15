// This file draws the velocity from Maxwell-Boltzmann distribution and evenly distribute the initial kinetic energy in each direction.
// CHECKED.

#include"Datastruct.h"
#include"functions.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

// Initialize the simulation field
// allocate memory for the simulation field and the coefficients.
// set the initial field with zeros.
// call the interpolation function for the coefficients.
// call the initiliazation of the velocity
int InitField(Sim *sim, Mat rho, Field *field, Coeff *coeff)
{
	int Nx, Ny, Nz;

	// retrieve the simulation domain size
	Nx = rho.Nx;
	Ny = rho.Ny;
	Nz = rho.Nz;
	field->Nx = Nx;
	field->Ny = Ny;
	field->Nz = Nz;

	// Set the sizes of the simulation field
	SetMatSize(&field->Sxx, Nx, Ny, Nz);
	SetMatSize(&field->Syy, Nx, Ny, Nz);
	SetMatSize(&field->Szz, Nx, Ny, Nz);
	SetMatSize(&field->Syz, Nx, Ny, Nz);
	SetMatSize(&field->Sxz, Nx, Ny, Nz);
	SetMatSize(&field->Sxy, Nx, Ny, Nz);
	SetMatSize(&field->Vx, Nx, Ny, Nz);
	SetMatSize(&field->Vy, Nx, Ny, Nz);
	SetMatSize(&field->Vz, Nx, Ny, Nz);
	SetMatSize(&field->Vpx, Nx, Ny, Nz);
	SetMatSize(&field->Vpy, Nx, Ny, Nz);
	SetMatSize(&field->Vpz, Nx, Ny, Nz);
	// Allocate memory and initialize with zeros.
	MallocMat(&field->Sxx, 0.0);
	MallocMat(&field->Syy, 0.0);
	MallocMat(&field->Szz, 0.0);
	MallocMat(&field->Syz, 0.0);
	MallocMat(&field->Sxz, 0.0);
	MallocMat(&field->Sxy, 0.0);
	MallocMat(&field->Vx, 0.0);
	MallocMat(&field->Vy, 0.0);
	MallocMat(&field->Vz, 0.0);
	MallocMat(&field->Vpx, 0.0);
	MallocMat(&field->Vpy, 0.0);
	MallocMat(&field->Vpz, 0.0);

	// Set the size of the coefficient matrices
	SetMatSize(&coeff->buox, Nx, Ny, Nz);
	SetMatSize(&coeff->buoy, Nx, Ny, Nz);
	SetMatSize(&coeff->buoz, Nx, Ny, Nz);
	SetMatSize(&coeff->Muyz, Nx, Ny, Nz);
	SetMatSize(&coeff->Muxz, Nx, Ny, Nz);
	SetMatSize(&coeff->Muxy, Nx, Ny, Nz);

	// Allocate Memory and initialize with zeros.
	MallocMat(&coeff->buox, 0.0);
	MallocMat(&coeff->buoy, 0.0);
	MallocMat(&coeff->buoz, 0.0);
	MallocMat(&coeff->Muyz, 0.0);
	MallocMat(&coeff->Muxz, 0.0);
	MallocMat(&coeff->Muxy, 0.0);

	// Initialize RNG seed
	InitRand((unsigned long long)time(NULL));
	// Call the interpolation function for the coefficient matrices.
	InterpolateFreeCoeff(rho, coeff);
	// Call the initialization function for the Velocity.
	InitVelocity(sim, rho, field, coeff);
	// Give physically correct value to the coefficients.
	PhysicalCoeff(*sim, coeff);
	// Initialize sampling data set, hard code the length now.
	InitSampling(sim, sim->No);
	// normal exit
	return 0;
}
// Initialize velocity:
// 1. Draw the value from normal distribution with Box-Muller method.
// 2. Enforce zero net momentum in each direction.
// 3. All scale to correspond to nkT/m for the summation of the square of the velocities.

int InitVelocity(Sim *sim, Mat rho, Field *field, Coeff *coeff)
{
	int i, j, k;
	int Nx, Ny, Nz;
	double cx, cy, cz;
	double m1x, m1y, m1z;			// 1st order moment: momenum
	double m2x, m2y, m2z;			// 2nd order moment: kinetic energy
	double m3x, m3y, m3z;			// 3rd order moment: flux

	m1x = 0.0;
	m1y = 0.0;
	m1z = 0.0;
	cx = 0;
	cy = 0;
	cz = 0;

	sim->Natom = 0;
	// retrieve the size of the matrix.
	Nx = field->Nx;
	Ny = field->Ny;
	Nz = field->Nz;

	// first assign random number for the velocities.
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
			{
				// calculate the number of non zero mass density
				if (rho.e[i][j][k] != 0.0)
					sim->Natom++;

				// X direction
				if (coeff->buox.e[i][j][k] != 0)
				{
					field->Vx.e[i][j][k] = RandG(0,1);
					m1x += field->Vx.e[i][j][k]/ coeff->buox.e[i][j][k];
					cx += 1.0 / coeff->buox.e[i][j][k];
				}

				// Y direction
				if (coeff->buoy.e[i][j][k] != 0)
				{
					field->Vy.e[i][j][k] = RandG(0,1);
					m1y += field->Vy.e[i][j][k] / coeff->buoy.e[i][j][k];
					cy += 1.0 / coeff->buoy.e[i][j][k];
				}

				// Z direction
				if (coeff->buoz.e[i][j][k] != 0)
				{
					field->Vz.e[i][j][k] = RandG(0,1);
					m1z += field->Vz.e[i][j][k] / coeff->buoz.e[i][j][k];
					cz += 1.0 / coeff->buoz.e[i][j][k];
				}
			}

	// average of the mmt
	m1x /= cx;
	m1y /= cy;
	m1z /= cz;
	printf("Biased Vector = (%e, %e, %e)\n", m1x, m1y, m1z);

	// kinetic energy in each direction.
	m2x = 0.0;
	m2y = 0.0;
	m2z = 0.0;

	// compare the counted cx, cy and cz, they should be the same.
	printf("(Nx, Ny, Nz) = (%.2f, %.2f, %.2f), Number of cells: %d\n", cx, cy, cz, sim->Natom);

	// achieve strict zero net momentum
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
			{
				if (coeff->buox.e[i][j][k] != 0)
				{
					field->Vx.e[i][j][k] -= m1x;
					m2x += field->Vx.e[i][j][k] * field->Vx.e[i][j][k] / coeff->buox.e[i][j][k];
				}
				if (coeff->buoy.e[i][j][k] != 0)
				{
					field->Vy.e[i][j][k] -= m1y;
					m2y += field->Vy.e[i][j][k] * field->Vy.e[i][j][k] / coeff->buoy.e[i][j][k];
				}
				if (coeff->buoz.e[i][j][k] != 0)
				{
					field->Vz.e[i][j][k] -= m1z;
					m2z += field->Vz.e[i][j][k] * field->Vz.e[i][j][k] / coeff->buoz.e[i][j][k];
				}
			}

	// evenly distribute the kinetic energy, scale the velocity to correspond to n (According to equipartition theorem, should be nkT/m);
	// scaling factor.
	m2x = sqrt(sim->Natom / m2x);
	m2y = sqrt(sim->Natom / m2y);
	m2z = sqrt(sim->Natom / m2z);

	printf("Scaling Vector = (%e, %e, %e)\n", m2x, m2y, m2z);

	m1x = 0.0;
	m1y = 0.0;
	m1z = 0.0;
	m3x = 0.0;
	m3y = 0.0;
	m3z = 0.0;
	// scaling
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
			{
				field->Vx.e[i][j][k] *= m2x;
				field->Vy.e[i][j][k] *= m2y;
				field->Vz.e[i][j][k] *= m2z;
				// recalculate the initial net momentum
				if (coeff->buox.e[i][j][k] != 0)
				{
					m1x += field->Vx.e[i][j][k] / coeff->buox.e[i][j][k];
					m3x += field->Vx.e[i][j][k] * field->Vx.e[i][j][k] * field->Vx.e[i][j][k] / coeff->buox.e[i][j][k];
				}
				if (coeff->buoy.e[i][j][k] != 0)
				{
					m1y += field->Vy.e[i][j][k] / coeff->buoy.e[i][j][k];
					m3y += field->Vy.e[i][j][k] * field->Vy.e[i][j][k] * field->Vy.e[i][j][k] / coeff->buoy.e[i][j][k];
				}
				if (coeff->buoz.e[i][j][k] != 0)
				{
					m1z += field->Vz.e[i][j][k] / coeff->buoz.e[i][j][k];
					m3z += field->Vz.e[i][j][k] * field->Vz.e[i][j][k] * field->Vz.e[i][j][k] / coeff->buoz.e[i][j][k];
				}
			}

	printf("Net Momentum at Initialization = (%e, %e, %e)\n", m1x, m1y, m1z);
	printf("Net Heat Current at Initialization = (%e, %e, %e)\n", m3x, m3y, m3z);

	// output to check the distribution
	FILE *fp;
	fp = fopen("Vx.log","a+");
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
			{
				if (coeff->buox.e[i][j][k] != 0)
					fprintf(fp,"%e\n", field->Vx.e[i][j][k]);
			}
	fclose(fp);
	// normal exit
	return 0;
}

// initialize the sampling space for the heat current and energy
int InitSampling(Sim *sim, int len)
{
	int i;

	// set the parameters
	sim->energy.counter = 0;
	sim->energy.len = len;
	sim->flux.counter = 0;
	sim->flux.len = len;
	sim->convf.counter = 0;
	sim->convf.len = len;
	sim->EiSumAve = 0.0;

	// allocate memory for the logger
	sim->energy.ke = (double*)malloc(len*sizeof(double));
	sim->energy.se = (double*)malloc(len*sizeof(double));
	sim->energy.te = (double*)malloc(len*sizeof(double));

	sim->flux.x = (double*)malloc(len*sizeof(double));
	sim->flux.y = (double*)malloc(len*sizeof(double));
	sim->flux.z = (double*)malloc(len*sizeof(double));

	sim->convf.x = (double*)malloc(len*sizeof(double));
	sim->convf.y = (double*)malloc(len*sizeof(double));
	sim->convf.z = (double*)malloc(len*sizeof(double));

	// initialize with zeros
	for (i = 0; i < len; i++)
	{
		sim->energy.se[i] = 0.0;
		sim->energy.ke[i] = 0.0;
		sim->energy.te[i] = 0.0;

		sim->flux.x[i] = 0.0;
		sim->flux.y[i] = 0.0;
		sim->flux.z[i] = 0.0;

		sim->convf.x[i] = 0.0;
		sim->convf.y[i] = 0.0;
		sim->convf.z[i] = 0.0;
	}

	// normal exit
	return 0;
}
