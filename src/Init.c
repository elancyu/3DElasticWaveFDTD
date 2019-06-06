// This file draws the velocity from Maxwell-Boltzmann distribution and evenly distribute the initial kinetic energy in each direction.
#include"Datastruct.h"
#include"functions.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

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

	// Set the size of the coefficient matrices
	SetMatSize(&coeff->buox, Nx, Ny, Nz);
	SetMatSize(&coeff->buoy, Nx, Ny, Nz);
	SetMatSize(&coeff->buoz, Nx, Ny, Nz);
	SetMatSize(&coeff->Muyz, Nx, Ny, Nz);
	SetMatSize(&coeff->Muxz, Nx, Ny, Nz);
	SetMatSize(&coeff->Muxy, Nx, Ny, Nz);
	SetMatSize(&coeff->Lam2mu, Nx, Ny, Nz);
	SetMatSize(&coeff->Lam, Nx, Ny, Nz);

	// Allocate Memory and initialize with zeros.
	MallocMat(&coeff->buox, 0.0);
	MallocMat(&coeff->buoy, 0.0);
	MallocMat(&coeff->buoz, 0.0);
	MallocMat(&coeff->Muyz, 0.0);
	MallocMat(&coeff->Muxz, 0.0);
	MallocMat(&coeff->Muxy, 0.0);
	MallocMat(&coeff->Lam2mu, 0.0);
	MallocMat(&coeff->Lam, 0.0);

	// Call the interpolation function for the coefficient matrices.
	InterpolateCoeffs(rho, coeff);
	// Call the initialization function for the Velocity.
	InitVelocity(sim, rho, field, coeff);
	// Give physically correct value to the coefficients.
	PhysicalCoeff(*sim, coeff);
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
	double mmx, mmy, mmz;
	double cx, cy, cz;
	double kx, ky, kz;

	mmx = 0.0;
	mmy = 0.0;
	mmz = 0.0;
	cx = 0;
	cy = 0;
	cz = 0;


	sim->Natom = 0;
	// no need to check the size of the matrix, assuming periodic boundary condition in X and Y directions.

	Nx = field->Nx;
	Ny = field->Ny;
	Nz = field->Nz;

	// count the number of nonzero mass density
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
				if (rho.e[i][j][k] != 0.0)
					sim->Natom++;

	// first assign random number for the velocities.
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
			{
				// X direction
				if (coeff->buox.e[i][j][k] != 0)
				{
					field->Vx.e[i][j][k] = RandG(0,1);
					mmx += field->Vx.e[i][j][k] / coeff->buox.e[i][j][k];
					cx += 1.0 / coeff->buox.e[i][j][k];
				}

				// Y direction
				if (coeff->buoy.e[i][j][k] != 0)
				{
					field->Vy.e[i][j][k] = RandG(0,1);
					mmy += field->Vy.e[i][j][k] / coeff->buoy.e[i][j][k];
					cy += 1.0 / coeff->buoy.e[i][j][k];
				}

				// Z direction
				if (coeff->buoz.e[i][j][k] != 0)
				{
					field->Vz.e[i][j][k] = RandG(0,1);
					mmz += field->Vz.e[i][j][k] / coeff->buoz.e[i][j][k];
					cz += 1.0 / coeff->buoz.e[i][j][k];
				}
			}

	// average of the mmt
	mmx /= cx;
	mmy /= cy;
	mmz /= cz;

	kx = 0.0;
	ky = 0.0;
	kz = 0.0;

	// compare the counted cx, cy and cz, they should be the same.
	printf("(Nx, Ny, Nz) = (%f, %f, %f), Number of cells: %d\n", cx, cy, cz, sim->Natom);

	// achieve strict zero net momentum
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
			{
				if (coeff->buox.e[i][j][k] != 0)
				{
					field->Vx.e[i][j][k] -= mmx / coeff->buox.e[i][j][k];
					kx += field->Vx.e[i][j][k] * field->Vx.e[i][j][k] / coeff->buox.e[i][j][k];
				}
				if (coeff->buoy.e[i][j][k] != 0)
				{
					field->Vy.e[i][j][k] -= mmy / coeff->buoy.e[i][j][k];
					ky += field->Vy.e[i][j][k] * field->Vy.e[i][j][k] / coeff->buoy.e[i][j][k];
				}
				if (coeff->buoz.e[i][j][k] != 0)
				{
						field->Vz.e[i][j][k] -= mmz / coeff->buoz.e[i][j][k];
						kz += field->Vz.e[i][j][k] / coeff->buoz.e[i][j][k];
				}
			}

	// evenly distribute the kinetic energy, scale the velocity to correspond to n (According to equipartition theorem, should be nkT/m);
	// scaling factor.
	kx = sim->Natom / kx;
	ky = sim->Natom / ky;
	kz = sim->Natom / kz;

	// scaling
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
			{
				field->Vx.e[i][j][k] *= kx;
				field->Vy.e[i][j][k] *= ky;
				field->Vz.e[i][j][k] *= kz;
			}

	// may need to output the initial velocity for the simulation?

	// normal exit
	return 0;
}


// Still need to take care of this part.
// interpolating the elastic constants for free boundary conditions
int InterpolateCoeffs(Mat rho, Coeff *coeff)
{
	FILE *fp;
	int i, j, k;
	int im, jm, km;
	int Nx, Ny, Nz;

	Nx = rho.Nx;
	Ny = rho.Ny;
	Nz = rho.Nz;

	// interpolate the elastic constants for mu
	// Mu_yz & Mu_xz
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 1; k < Nz; k++)
			{
				im = (i - 1 + Nx) % Nx;
				jm = (j - 1 + Ny) % Ny;
				km = k - 1;
				if (rho.e[i][j][k] + rho.e[i][jm][k] + rho.e[i][j][km] + rho.e[i][jm][km] == 4.0)
					coeff->Muyz.e[i][j][k] = 1.0;
				if (rho.e[i][j][k] + rho.e[im][j][k] + rho.e[i][j][km] + rho.e[im][j][km] == 4.0)
					coeff->Muxz.e[i][j][k] = 1.0;
			}

	// mu_xy
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
			{
				im = (i - 1 + Nx) % Nx;
				jm = (j - 1 + Ny) % Ny;
				if (rho.e[i][j][k] + rho.e[im][j][k] + rho.e[i][jm][k] + rho.e[im][jm][k] == 4.0)
					coeff->Muxy.e[i][j][k] = 1.0;
			}

	// interpolate elastic constant lam2mu
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
				if (rho.e[i][j][k] == 1.0)
				{
					coeff->Lam.e[i][j][k] = 1.0;
					coeff->Lam2mu.e[i][j][k] = 1.0;
				}

	// output the interpolated values
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
			{
				im = (i - 1 + Nx) % Nx;
				jm = (j - 1 + Ny) % Ny;
				km = (k - 1 + Nz) % Nz;
				// needs modification since periodic boundary condition is used in X and Y directions.
				if (rho.e[i][j][k] + rho.e[im][j][k] > 0.0)
					coeff->buox.e[i][j][k] = 2.0 / (rho.e[i][j][k] + rho.e[im][j][k]);
				if (rho.e[i][j][k] + rho.e[i][jm][k] > 0.0)
					coeff->buoy.e[i][j][k] = 2.0 / (rho.e[i][j][k] + rho.e[i][jm][k]);
				if (rho.e[i][j][k] + rho.e[i][j][km] > 0.0)
					coeff->buoz.e[i][j][k] = 2.0 / (rho.e[i][j][k] + rho.e[i][j][km]);
			}
	// record the interpolated coefficients
	fp = fopen("InterpolatedElasticConstants.log","a+");
	for (k = 0; k < Nz; k++)
		for (i = 0; i < Nx; i++)
			for (j = 0; j < Ny; j++)
			{
				fprintf(fp,"%e %e %e %e %e\n", coeff->Muyz.e[i][j][k], coeff->Muxz.e[i][j][k], coeff->Muxy.e[i][j][k], coeff->Lam.e[i][j][k], coeff->Lam2mu.e[i][j][k]);
			}
	fclose(fp);

	fp = fopen("InterpolatedBuoyancy.log","a+");
	for (k = 0; k < Nz; k++)
		for (i = 0; i < Nx; i++)
			for (j = 0; j < Ny; j++)
			{
				fprintf(fp,"%e %e %e\n", coeff->buox.e[i][j][k], coeff->buoy.e[i][j][k], coeff->buoz.e[i][j][k]);
			}
	fclose(fp);
	// normal exit
	return 0;
}

int PhysicalCoeff(Sim sim, Coeff *coeff)
{
	int i, j, k;
	int Nx, Ny, Nz;

	Nx = coeff->buox.Nx;
	Ny = coeff->buox.Ny;
	Nz = coeff->buox.Nz;

	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
			{
				coeff->buox.e[i][j][k] *= 1.0 / sim.rho;
				coeff->buoy.e[i][j][k] *= 1.0 / sim.rho;
				coeff->buoz.e[i][j][k] *= 1.0 / sim.rho;
				coeff->Muxy.e[i][j][k] *= sim.mu;
				coeff->Muxz.e[i][j][k] *= sim.mu;
				coeff->Muyz.e[i][j][k] *= sim.mu;
				coeff->Lam.e[i][j][k] *= sim.lambda;
				coeff->Lam2mu.e[i][j][k] *= (sim.lambda + 2*sim.mu);
			}
	// normal exit
	return 0;
}
