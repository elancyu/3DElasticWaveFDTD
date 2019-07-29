// This file contains the energy, momentum and flux calculations.
// Checked. The code is safe.

#include"functions.h"
#include"Datastruct.h"
#include<stdio.h>

// Strictly calculate the kinetic energy and strain energy.
int SampleEnergy(Sim *sim, Mat rho,Field field)
{
	int i, j, k;
	int im, jm, km;
	int ip, jp, kp;
	double wyz, wxz, wxy;
	double wx, wy, wz;
	double Ek, Es;
	double ek, es, e;
	int Nx, Ny, Nz;
	double vx, vy, vz;
	double evx, evy, evz;
	// remember to remove s11, s12, s44;
	Es = 0.0;
	Ek = 0.0;
	vx = 0.0;
	vy = 0.0;
	vz = 0.0;
	evx = 0.0;
	evy = 0.0;
	evz = 0.0;

	// retrieve the matrix size
	Nx = rho.Nx;
	Ny = rho.Ny;
	Nz = rho.Nz;

	// directly add up all Sii^2, and the cross term, SiiSjj.
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
			{
				if (rho.e[i][j][k] == 0)
					continue;

				// minus direction indices
				im = (i - 1 + Nx) % Nx;
				jm = (j - 1 + Ny) % Ny;
				km = (k - 1 + Nz) % Nz;
				ip = (i + 1) % Nx;
				jp = (j + 1) % Ny;
				kp = (k + 1) % Nz;
				es = 0.0;
				ek = 0.0;
				// Sii^2 && SiiSjj
				es = 0.5 * sim->s11 * field.Sxx.e[i][j][k] * field.Sxx.e[i][j][k];
				es += 0.5 * sim->s11 * field.Syy.e[i][j][k] * field.Syy.e[i][j][k];
				es += 0.5 * sim->s11 * field.Szz.e[i][j][k] * field.Szz.e[i][j][k];
				es += sim->s12 * field.Sxx.e[i][j][k] * field.Syy.e[i][j][k];
				es += sim->s12 * field.Sxx.e[i][j][k] * field.Szz.e[i][j][k];
				es += sim->s12 * field.Syy.e[i][j][k] * field.Szz.e[i][j][k];

				// Sij^2
				es += 0.25 * sim->s44 * (field.Syz.e[i][j][k] * field.Syz.e[i][j][k] + field.Syz.e[i][jp][k] * field.Syz.e[i][jp][k] + field.Syz.e[i][j][kp] * field.Syz.e[i][j][kp] + field.Syz.e[i][jp][kp] * field.Syz.e[i][jp][kp]);
				es += 0.25 * sim->s44 * (field.Sxz.e[i][j][k] * field.Sxz.e[i][j][k] + field.Sxz.e[ip][j][k] * field.Sxz.e[ip][j][k] + field.Sxz.e[i][j][kp] * field.Sxz.e[i][j][kp] + field.Sxz.e[ip][j][kp] * field.Sxz.e[ip][j][kp]);
				es += 0.25 * sim->s44 * (field.Sxy.e[i][j][k] * field.Sxy.e[i][j][k] + field.Sxy.e[ip][j][k] * field.Sxy.e[ip][j][k] + field.Sxy.e[i][jp][k] * field.Sxy.e[i][jp][k] + field.Sxy.e[ip][jp][k] * field.Sxy.e[ip][jp][k]);

				// Kinetic Energy
				ek = 0.25 * sim->rho * (field.Vx.e[i][j][k] * field.Vpx.e[i][j][k] + field.Vx.e[ip][j][k] * field.Vpx.e[ip][j][k]);
				ek += 0.25 * sim->rho * (field.Vy.e[i][j][k] * field.Vpy.e[i][j][k] + field.Vy.e[i][jp][k] * field.Vpy.e[i][jp][k]);
				ek += 0.25 * sim->rho * (field.Vz.e[i][j][k] * field.Vpz.e[i][j][k] + field.Vz.e[i][j][kp] * field.Vpz.e[i][j][kp]);

				// Ek and Es summation
				Es += es;
				Ek += ek;

				// extra information
				vx += 0.25 * (field.Vx.e[i][j][k] + field.Vx.e[ip][j][k] + field.Vpx.e[i][j][k] + field.Vpx.e[ip][j][k]);
				vy += 0.25 * (field.Vy.e[i][j][k] + field.Vy.e[i][jp][k] + field.Vpy.e[i][j][k] + field.Vpy.e[i][jp][k]);
				vz += 0.25 * (field.Vz.e[i][j][k] + field.Vz.e[i][j][kp] + field.Vpz.e[i][j][k] + field.Vpz.e[i][j][kp]);

				// Sigma eivi
				e = es + ek;
				evx += 0.25 * e * (field.Vx.e[i][j][k] + field.Vx.e[ip][j][k] + field.Vpx.e[i][j][k] + field.Vpx.e[ip][j][k]);
				evy += 0.25 * e * (field.Vy.e[i][j][k] + field.Vy.e[i][jp][k] + field.Vpy.e[i][j][k] + field.Vpy.e[i][jp][k]);
				evz += 0.25 * e * (field.Vz.e[i][j][k] + field.Vz.e[i][j][kp] + field.Vpz.e[i][j][k] + field.Vpz.e[i][j][kp]);
			}

	// times the volume
	Es *= sim->dV;
	Ek *= sim->dV;

	i = sim->energy.counter;
	sim->energy.ke[i] = Ek;
	sim->energy.se[i] = Es;
	sim->energy.te[i] = Ek + Es;

	i = sim->convf.counter;
	//printf("VSum = (%e, %e, %e)\n", vx, vy, vz);
	sim->convf.x[i] = vx;
	sim->convf.y[i] = vy;
	sim->convf.z[i] = vz;
	sim->convf.counter++;

	// extra information
	sim->evx = evx * sim->dV;
	sim->evy = evy * sim->dV;
	sim->evz = evz * sim->dV;

	// increment of the counter
	sim->energy.counter++;

	// if it is full, output and reset.
	if (sim->energy.counter == sim->energy.len)
		WriteEnergy(sim);
	// normal exit
	return 0;
}

// Sample the heat current spatial integral, first interpolate to velocity and then integrate over the domain
int SampleHeatCurrent(Sim *sim, Mat rho, Field field)
{
	int i, j, k;
	int im, jm, km, ip, jp, kp;			// interpolating indices;
	int Nx, Ny, Nz;
	double Jx, Jy, Jz;
	double VxSxx, VySxy, VzSxz;
	double VxSxy, VySyy, VzSyz;
	double VxSxz, VySyz, VzSzz;
	double Vx, Vy, Vz;
	double Sxx, Sxy, Sxz, Syy, Szz, Syz;
	double vx, vy, vz;

	Nx = rho.Nx;
	Ny = rho.Ny;
	Nz = rho.Nz;

	Jx = 0.0;
	Jy = 0.0;
	Jz = 0.0;

	vx = 0.0;
	vy = 0.0;
	vz = 0.0;

	/*
	// M2:integration method: reasoned directly by differencing the energy: seem to be fine.
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
			{
				// skip the virtual cells
				if (rho.e[i][j][k] == 0)
					continue;

				im = (i - 1 + Nx) % Nx;
				jm = (j - 1 + Ny) % Ny;
				km = (k - 1 + Nz) % Nz;
				ip = (i + 1) % Nx;
				jp = (j + 1) % Ny;
				kp = (k + 1) % Nz;

				// x direction
				// VxSxx: interpolate to velocity (X- surface)
				VxSxx = 0.25 * (field.Vx.e[i][j][k] + field.Vpx.e[i][j][k]) * (field.Sxx.e[i][j][k] + field.Sxx.e[im][j][k]);

				// VySxy: interpolate to the stress location
				VySxy = 0.125 * field.Sxy.e[i][j][k] * (field.Vy.e[im][j][k] + field.Vy.e[i][j][k] + field.Vpy.e[im][j][k] + field.Vpy.e[i][j][k]);
				VySxy += 0.125 * field.Sxy.e[i][jp][k] * (field.Vy.e[im][jp][k] + field.Vy.e[i][jp][k] + field.Vpy.e[im][jp][k] + field.Vpy.e[i][jp][k]);

				// VzSxz: inteprolate to the stress location
				VzSxz = 0.125 * field.Sxz.e[i][j][k] * (field.Vz.e[im][j][k] + field.Vz.e[i][j][k] + field.Vpz.e[im][j][k] + field.Vpz.e[i][j][k]);
				VzSxz += 0.125 * field.Sxz.e[i][j][kp] * (field.Vz.e[im][j][kp] + field.Vz.e[i][j][kp] + field.Vpz.e[im][j][kp] + field.Vpz.e[i][j][kp]);
				// Jx
				Jx -= (VxSxx + VySxy + VzSxz) * rho.e[im][j][k];

				// y direction
				// VxSxy
				VxSxy = 0.125 * field.Sxy.e[i][j][k] * (field.Vx.e[i][jm][k] + field.Vx.e[i][j][k] + field.Vpx.e[i][jm][k] + field.Vpx.e[i][j][k]);
				VxSxy += 0.125 * field.Sxy.e[ip][j][k] * (field.Vx.e[ip][jm][k] + field.Vx.e[ip][j][k] + field.Vpx.e[ip][jm][k] + field.Vpx.e[ip][j][k]);
				// VySyy
				VySyy = 0.25 * (field.Vy.e[i][j][k] + field.Vpy.e[i][j][k]) * (field.Syy.e[i][j][k] + field.Syy.e[i][jm][k]);

				// VzSyz
				VzSyz = 0.125 * field.Syz.e[i][j][k] * (field.Vz.e[i][jm][k] + field.Vz.e[i][j][k] + field.Vpz.e[i][jm][k] + field.Vpz.e[i][j][k]);
				VzSyz += 0.125 * field.Syz.e[i][j][kp] * (field.Vz.e[i][jm][kp] + field.Vz.e[i][j][kp] + field.Vpz.e[i][jm][kp] + field.Vpz.e[i][j][kp]);
				// jy
				Jy -= (VxSxy + VySyy + VzSyz) * rho.e[i][jm][k];

				// z direction
				// VxSxz
				VxSxz = 0.125 * field.Sxz.e[i][j][k] * (field.Vx.e[i][j][km] + field.Vx.e[i][j][k] + field.Vpx.e[i][j][km] + field.Vpx.e[i][j][k]);
				VxSxz += 0.125 * field.Sxz.e[ip][j][k] * (field.Vx.e[ip][j][km] + field.Vx.e[ip][j][k] + field.Vpx.e[ip][j][km] + field.Vpx.e[ip][j][k]);
				// VySyz
				VySyz = 0.125 * field.Syz.e[i][j][k] * (field.Vy.e[i][j][km] + field.Vy.e[i][j][k] + field.Vpy.e[i][j][km] + field.Vpy.e[i][j][k]);
				VySyz += 0.25 * field.Syz.e[i][jp][k] * (field.Vy.e[i][jp][km] + field.Vy.e[i][jp][k] + field.Vpy.e[i][jp][km] + field.Vpy.e[i][jp][k]);
				// VzSzz
				VzSzz = 0.25 * (field.Vz.e[i][j][k] + field.Vpz.e[i][j][k]) * (field.Szz.e[i][j][k] + field.Szz.e[i][j][km]);
				// jz
				Jz -= (VxSxz + VySyz + VzSzz) * rho.e[i][j][km];
			}
			*/

	// M2: first interpolate the physical quantities to the surface center and then do the integration.
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
			{
				if (rho.e[i][j][k] == 0)
					continue;

				im = (i - 1 + Nx) % Nx;
				jm = (j - 1 + Ny) % Ny;
				km = (k - 1 + Nz) % Nz;
				ip = (i + 1) % Nx;
				jp = (j + 1) % Ny;
				kp = (k + 1) % Nz;

				// interpolation to X- surface.
				Vx = 0.5 * (field.Vx.e[i][j][k] + field.Vpx.e[i][j][k]);
				Vy = 0.125 * (field.Vy.e[i][j][k] + field.Vy.e[im][j][k] + field.Vy.e[i][jp][k] + field.Vy.e[im][jp][k] + field.Vpy.e[i][j][k] + field.Vpy.e[im][j][k] + field.Vpy.e[i][jp][k] + field.Vpy.e[im][jp][k]);
				Vz = 0.125 * (field.Vz.e[i][j][k] + field.Vz.e[im][j][k] + field.Vz.e[i][j][kp] + field.Vz.e[im][j][kp] + field.Vpz.e[i][j][k] + field.Vpz.e[im][j][k] + field.Vpz.e[i][j][kp] + field.Vpz.e[im][j][kp]);
				Sxx = 0.5 * (field.Sxx.e[i][j][k] + field.Sxx.e[im][j][k]);
				Sxy = 0.5 * (field.Sxy.e[i][j][k] + field.Sxy.e[i][jp][k]);
				Sxz = 0.5 * (field.Sxz.e[i][j][k] + field.Sxz.e[i][j][kp]);
				Jx -= (Vx * Sxx + Vy * Sxy + Vz * Sxz) * rho.e[im][j][k];			// if it is at the interface, remove it.

				// interpolation to Y- surface
				Vx = 0.125 * (field.Vx.e[i][j][k] + field.Vx.e[i][jm][k] + field.Vx.e[ip][j][k] + field.Vx.e[ip][jm][k] + field.Vpx.e[i][j][k] + field.Vpx.e[i][jm][k] + field.Vpx.e[ip][j][k] + field.Vpx.e[ip][jm][k]);
				Vy = 0.5 * (field.Vy.e[i][j][k] + field.Vpy.e[i][j][k]);
				Vz = 0.125 * (field.Vz.e[i][j][k] + field.Vz.e[i][jm][k] + field.Vz.e[i][j][kp] + field.Vz.e[i][jm][kp] + field.Vpz.e[i][j][k] + field.Vpz.e[i][jm][k] + field.Vpz.e[i][j][kp] + field.Vpz.e[i][jm][kp]);
				Sxy = 0.5 * (field.Sxy.e[i][j][k] + field.Sxy.e[ip][j][k]);
				Syy = 0.5 * (field.Syy.e[i][j][k] + field.Syy.e[i][jm][k]);
				Syz = 0.5 * (field.Syz.e[i][j][k] + field.Syz.e[i][j][kp]);
				Jy -= (Vx * Sxy + Vy * Syy + Vz * Syz) * rho.e[i][jm][k];			//if it is at the surface, remove its contribution.

				// interpolation to Z- surface.
				Vx = 0.125 * (field.Vx.e[i][j][k] + field.Vx.e[i][j][km] + field.Vx.e[ip][j][k] + field.Vx.e[ip][j][km] + field.Vpx.e[i][j][k] + field.Vpx.e[i][j][km] + field.Vpx.e[ip][j][k] + field.Vpx.e[ip][j][km]);
				Vy = 0.125 * (field.Vy.e[i][j][k] + field.Vy.e[i][j][km] + field.Vy.e[i][jp][k] + field.Vy.e[i][jp][km] + field.Vpy.e[i][j][k] + field.Vpy.e[i][j][km] + field.Vpy.e[i][jp][k] + field.Vpy.e[i][jp][km]);
				Vz = 0.5 * (field.Vz.e[i][j][k] + field.Vpz.e[i][j][k]);
				Sxz = 0.5 * (field.Sxz.e[i][j][k] + field.Sxz.e[ip][j][k]);
				Syz = 0.5 * (field.Syz.e[i][j][k] + field.Syz.e[i][jp][k]);
				Szz = 0.5 * (field.Szz.e[i][j][k] + field.Szz.e[i][j][km]);
				Jz -= (Vx * Sxz + Vy * Syz + Vz * Szz) * rho.e[i][j][km];			// if it is at the interface, remove it.
			}

	Jx *= sim->dV;
	Jy *= sim->dV;
	Jz *= sim->dV;
	// 5th method: Gaussian Quadrature, first interpolate to Gaussian point, then multiply and

	i = sim->flux.counter;
	// sim->flux.x[i] = sim->evx - Jx;
	// sim->flux.y[i] = sim->evy - Jy;
	// sim->flux.z[i] = sim->evz - Jz;
	sim->flux.x[i] = Jx;
	sim->flux.y[i] = Jy;
	sim->flux.z[i] = Jz;
	sim->flux.counter++;
	// if it is full, output and reset.
	if (sim->flux.counter == sim->flux.len)
		WriteHeatCurrent(sim);
	// normal exit
	return 0;
}

// output the heat current
int WriteHeatCurrent(Sim *sim)
{
	FILE *fp;
	int i, N;
	N = sim->flux.len;

	fp = fopen("HC.out","a+");
	for (i = 0; i < N; i++)
	{
		fprintf(fp, "%e %e %e\n",sim->flux.x[i], sim->flux.y[i], sim->flux.z[i]);
		// reset to zeros
		sim->flux.x[i] = 0.0;
		sim->flux.y[i] = 0.0;
		sim->flux.z[i] = 0.0;
	}
	// reset the counter
	sim->flux.counter = 0;
	fclose(fp);


	// normal exit
	return 0;
}

// Set the Heat Current and the Energy Records to zeros
int WriteEnergy(Sim *sim)
{
	FILE *fp;
	int i;
	int N;
	N = sim->energy.len;

	fp = fopen("Energy.log","a+");
	for (i = 0; i < N; i++)
	{
		fprintf(fp,"%e %e %e\n", sim->energy.ke[i], sim->energy.se[i], sim->energy.te[i]);
		sim->energy.ke[i] = 0.0;
		sim->energy.se[i] = 0.0;
		sim->energy.te[i] = 0.0;
	}
	sim->energy.counter = 0;
	fclose(fp);

	fp = fopen("VSum.out","a+");
	for (i = 0; i < N; i++)
	{
		fprintf(fp,"%.5e %.5e %.5e\n", sim->convf.x[i], sim->convf.y[i], sim->convf.z[i]);
		sim->convf.x[i] = 0.0;
		sim->convf.y[i] = 0.0;
		sim->convf.z[i] = 0.0;
	}
	sim->convf.counter = 0;
	fclose(fp);
	// normal exit
	return 0;
}
