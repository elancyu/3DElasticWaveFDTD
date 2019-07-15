// This file contains the coefficient preparation functions.
#include"Datastruct.h"
#include"functions.h"
#include<stdio.h>
#include<stdlib.h>

// Still need to take care of this part.
// interpolating the elastic constants for free boundary conditions
int InterpolateFreeCoeff(Mat rho, Coeff *coeff)
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
			for (k = 0; k < Nz; k++)
			{
				im = (i - 1 + Nx) % Nx;
				jm = (j - 1 + Ny) % Ny;
				km = (k - 1 + Nz) % Nz;
				// elastic constant Mu
				if (rho.e[i][j][k] + rho.e[i][jm][k] + rho.e[i][j][km] + rho.e[i][jm][km] == 4.0)
					coeff->Muyz.e[i][j][k] = 1.0;
				if (rho.e[i][j][k] + rho.e[im][j][k] + rho.e[i][j][km] + rho.e[im][j][km] == 4.0)
					coeff->Muxz.e[i][j][k] = 1.0;
				if (rho.e[i][j][k] + rho.e[im][j][k] + rho.e[i][jm][k] + rho.e[im][jm][k] == 4.0)
					coeff->Muxy.e[i][j][k] = 1.0;

				// buoyancies
				// needs modification since periodic boundary condition is used in X and Y directions.
				if (rho.e[i][j][k] + rho.e[im][j][k] > 0.0)
					coeff->buox.e[i][j][k] = 2.0 / (rho.e[i][j][k] + rho.e[im][j][k]);
				if (rho.e[i][j][k] + rho.e[i][jm][k] > 0.0)
					coeff->buoy.e[i][j][k] = 2.0 / (rho.e[i][j][k] + rho.e[i][jm][k]);
				if (rho.e[i][j][k] + rho.e[i][j][km] > 0.0)
					coeff->buoz.e[i][j][k] = 2.0 / (rho.e[i][j][k] + rho.e[i][j][km]);
			}

	// record the interpolated coefficients
	fp = fopen("FreeBCElasticConstants.log","a+");
	for (k = 0; k < Nz; k++)
		for (i = 0; i < Nx; i++)
			for (j = 0; j < Ny; j++)
				fprintf(fp,"%e %e %e\n", coeff->Muyz.e[i][j][k], coeff->Muxz.e[i][j][k], coeff->Muxy.e[i][j][k]);
	fclose(fp);

	fp = fopen("FreeBCBuoyancy.log","a+");
	for (k = 0; k < Nz; k++)
		for (i = 0; i < Nx; i++)
			for (j = 0; j < Ny; j++)
				fprintf(fp,"%e %e %e\n", coeff->buox.e[i][j][k], coeff->buoy.e[i][j][k], coeff->buoz.e[i][j][k]);
	fclose(fp);
	// normal exit
	return 0;
}

// Still need to take care of this part.
// interpolating the elastic constants for fixed boundary conditions
int InterpolateFixedCoeff(Mat rho, Coeff *coeff)
{
	FILE *fp;
	int i, j, k;
	int im, jm, km;
	int Nx, Ny, Nz;
	double rhosum;

	Nx = rho.Nx;
	Ny = rho.Ny;
	Nz = rho.Nz;

	// interpolate the elastic constants for mu
	// Mu_yz & Mu_xz
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
			{
				im = (i - 1 + Nx) % Nx;
				jm = (j - 1 + Ny) % Ny;
				km = (k - 1 + Nz) % Nz;
				// Mu_yz
				rhosum = rho.e[i][j][k] + rho.e[i][jm][k] + rho.e[i][j][km] + rho.e[i][jm][km];
				if (rhosum == 4.0 || rhosum == 3.0)
					coeff->Muyz.e[i][j][k] = 1.0;
				else if (rhosum == 2.0)
					coeff->Muyz.e[i][j][k] = 2.0;
				// Mu_xz
				rhosum = rho.e[i][j][k] + rho.e[im][j][k] + rho.e[i][j][km] + rho.e[im][j][km];
				if (rhosum == 4.0 || rhosum == 3.0)
					coeff->Muxz.e[i][j][k] = 1.0;
				else if (rhosum == 2.0)
					coeff->Muxz.e[i][j][k] = 2.0;
				// Mu_xy
				rhosum = rho.e[i][j][k] + rho.e[im][j][k] + rho.e[i][jm][k] + rho.e[im][jm][k];
				if (rhosum == 4.0 || rhosum == 3.0)
					coeff->Muxy.e[i][j][k] = 1.0;
				else if (rhosum == 2.0)
					coeff->Muxy.e[i][j][k] = 2.0;
				// buoyancies
				if (rho.e[i][j][k] * rho.e[im][j][k] > 0.0)
					coeff->buox.e[i][j][k] = 2.0 / (rho.e[i][j][k] + rho.e[im][j][k]);
				if (rho.e[i][j][k] * rho.e[i][jm][k] > 0.0)
					coeff->buoy.e[i][j][k] = 2.0 / (rho.e[i][j][k] + rho.e[i][jm][k]);
				if (rho.e[i][j][k] * rho.e[i][j][km] > 0.0)
					coeff->buoz.e[i][j][k] = 2.0 / (rho.e[i][j][k] + rho.e[i][j][km]);
			}

	// record the interpolated coefficients
	fp = fopen("FixedBCElasticConstants.log","a+");
	for (k = 0; k < Nz; k++)
		for (i = 0; i < Nx; i++)
			for (j = 0; j < Ny; j++)
				fprintf(fp,"%e %e %e\n", coeff->Muyz.e[i][j][k], coeff->Muxz.e[i][j][k], coeff->Muxy.e[i][j][k]);
	fclose(fp);

	fp = fopen("FixedBCBuoyancy.log","a+");
	for (k = 0; k < Nz; k++)
		for (i = 0; i < Nx; i++)
			for (j = 0; j < Ny; j++)
				fprintf(fp,"%e %e %e\n", coeff->buox.e[i][j][k], coeff->buoy.e[i][j][k], coeff->buoz.e[i][j][k]);
	fclose(fp);
	// normal exit
	return 0;
}

int PhysicalCoeff(Sim sim, Coeff *coeff)
{
	int i, j, k;
	int Nx, Ny, Nz;
	double buo;

	buo = 1.0 / sim.rho;

	Nx = coeff->buox.Nx;
	Ny = coeff->buox.Ny;
	Nz = coeff->buox.Nz;

	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
			{
				coeff->buox.e[i][j][k] *= buo;
				coeff->buoy.e[i][j][k] *= buo;
				coeff->buoz.e[i][j][k] *= buo;
				coeff->Muxy.e[i][j][k] *= sim.mu;
				coeff->Muxz.e[i][j][k] *= sim.mu;
				coeff->Muyz.e[i][j][k] *= sim.mu;
			}
	// normal exit
	return 0;
}
