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
	// remember to remove s11, s12, s44;
	Es = 0.0;
	Ek = 0.0;

	// retrieve the matrix size
	Nx = rho.Nx;
	Ny = rho.Ny;
	Nz = rho.Nz;

	// eivi
	sim->evx = 0.0;
	sim->evy = 0.0;
	sim->evz = 0.0;
	// directly add up all Sii^2, and the cross term, SiiSjj.
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
			{
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
				es = 0.5 * sim->s11 * rho.e[i][j][k] * field.Sxx.e[i][j][k] * field.Sxx.e[i][j][k];
				es += 0.5 * sim->s11 * rho.e[i][j][k] * field.Syy.e[i][j][k] * field.Syy.e[i][j][k];
				es += 0.5 * sim->s11 * rho.e[i][j][k] * field.Szz.e[i][j][k] * field.Szz.e[i][j][k];
				es += sim->s12 * rho.e[i][j][k] * field.Sxx.e[i][j][k] * field.Syy.e[i][j][k];
				es += sim->s12 * rho.e[i][j][k] * field.Sxx.e[i][j][k] * field.Szz.e[i][j][k];
				es += sim->s12 * rho.e[i][j][k] * field.Syy.e[i][j][k] * field.Szz.e[i][j][k];
				// weighting factor for Sij^2
				wyz = 0.25 * (rho.e[i][j][k] + rho.e[i][jm][k] + rho.e[i][j][km] + rho.e[i][jm][km]);
				wxz = 0.25 * (rho.e[i][j][k] + rho.e[im][j][k] + rho.e[i][j][km] + rho.e[im][j][km]);
				wxy = 0.25 * (rho.e[i][j][k] + rho.e[im][j][k] + rho.e[i][jm][k] + rho.e[im][jm][k]);
				// Sij^2
				es += sim->s44 * wyz * field.Syz.e[i][j][k] * field.Syz.e[i][j][k];
				es += sim->s44 * wxz * field.Sxz.e[i][j][k] * field.Sxz.e[i][j][k];
				es += sim->s44 * wxy * field.Sxy.e[i][j][k] * field.Sxy.e[i][j][k];

				// Kinetic Energy
				wx = 0.5 * (rho.e[i][j][k] + rho.e[im][j][k]);
				wy = 0.5 * (rho.e[i][j][k] + rho.e[i][jm][k]);
				wz = 0.5 * (rho.e[i][j][k] + rho.e[i][j][km]);
				ek = 0.5 * wx * sim->rho * field.Vx.e[i][j][k] * field.Vx.e[i][j][k];
				ek += 0.5 * wy * sim->rho * field.Vy.e[i][j][k] * field.Vy.e[i][j][k];
				ek += 0.5 *wz * sim->rho * field.Vz.e[i][j][k] * field.Vz.e[i][j][k];

				// Ek and Es summation
				Es += es;
				Ek += ek;

				// eivi
				e = es + ek;
				sim->evx += 0.5 * e * (field.Vx.e[i][j][k] + field.Vx.e[ip][j][k]);
				sim->evy += 0.5 * e * (field.Vy.e[i][j][k] + field.Vy.e[i][jp][k]);
				sim->evz += 0.5 * e * (field.Vz.e[i][j][k] + field.Vz.e[i][j][kp]);
			}

	// times the volume
	Es *= sim->dV;
	Ek *= sim->dV;

	i = sim->energy.counter;
	sim->energy.ke[i] = Ek;
	sim->energy.se[i] = Es;
	sim->energy.te[i] = Ek + Es;
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
	double VxSxx, VySxy, VzSxz; 		// components of Jx
	double VxSxy, VySyy, VzSyz;			// components of Jy
	double VxSxz, VySyz, VzSzz;			// components of Jz
	double jx, jy, jz;
	double sxx, syy, szz;
	double syz, sxz, sxy;

	Nx = rho.Nx;
	Ny = rho.Ny;
	Nz = rho.Nz;

	jx = 0.0;
	jy = 0.0;
	jz = 0.0;

	/*// 1st integration method: first interpolate the stress to the location of the velocity and
	// then calculate the heat current component and then interpolate back into the cell center.
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

				// VxSxx
				VxSxx = field.Vx.e[i][j][k] * (field.Sxx.e[i][j][k] + field.Sxx.e[im][j][k]);
				VxSxx += field.Vx.e[ip][j][k] * (field.Sxx.e[i][j][k] + field.Sxx.e[ip][j][k]);
				VxSxx /= -4;
				// VySxy
				VySxy = field.Vy.e[i][j][k] * (field.Sxy.e[i][j][k] + field.Sxy.e[ip][j][k]);
				VySxy += field.Vy.e[i][jp][k] * (field.Sxy.e[i][jp][k] + field.Sxy.e[ip][jp][k]);
				VySxy /= -4;
				// VzSxz
				VzSxz = field.Vz.e[i][j][k] * (field.Sxz.e[i][j][k] + field.Sxz.e[ip][j][k]);
				VzSxz += field.Vz.e[i][j][kp] * (field.Sxz.e[i][j][kp] + field.Sxz.e[ip][j][kp]);
				VzSxz /= -4;
				// Jx
				jx += (VxSxx + VySxy + VzSxz);

				// VxSxy
				VxSxy = field.Vx.e[i][j][k] * (field.Sxy.e[i][j][k] + field.Sxy.e[i][jp][k]);
				VxSxy += field.Vx.e[ip][j][k] * (field.Sxy.e[ip][j][k] + field.Sxy.e[ip][jp][k]);
				VxSxy /= -4;
				// VySyy
				VySyy = field.Vy.e[i][j][k] * (field.Syy.e[i][j][k] + field.Syy.e[i][jm][k]);
				VySyy += field.Vy.e[i][jp][k] * (field.Syy.e[i][j][k] + field.Syy.e[i][jp][k]);
				VySyy /= -4;
				// VzSyz
				VzSyz = field.Vz.e[i][j][k] * (field.Syz.e[i][j][k] + field.Syz.e[i][jp][k]);
				VzSyz += field.Vz.e[i][j][kp] * (field.Syz.e[i][j][kp] + field.Syz.e[i][jp][kp]);
				VzSyz /= -4;
				// Jy
				jy += (VxSxy + VySyy + VzSyz);

				// VxSxz
				VxSxz = field.Vx.e[i][j][k] * (field.Sxz.e[i][j][k] + field.Sxz.e[i][j][kp]);
				VxSxz += field.Vx.e[ip][j][k] * (field.Sxz.e[ip][j][k] + field.Sxz.e[ip][j][kp]);
				VxSxz /= -4;
				// VySyz
				VySyz = field.Vy.e[i][j][k] * (field.Syz.e[i][j][k] + field.Syz.e[i][j][kp]);
				VySyz += field.Vy.e[i][jp][k] * (field.Syz.e[i][jp][k] + field.Syz.e[i][jp][kp]);
				VySyz /= -4;
				// VzSzz
				VzSzz = field.Vz.e[i][j][k] * (field.Szz.e[i][j][k] + field.Szz.e[i][j][km]);
				VzSzz += field.Vz.e[i][j][kp] * (field.Szz.e[i][j][k] + field.Szz.e[i][j][kp]);
				VzSzz /= -4;
				// Jz
				jz += (VxSxz + VySyz + VzSzz);
			} */

	// 2nd integration method: reasoned directly by differencing the energy.
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
				// VxSxx
				VxSxx = (field.Vx.e[i][j][k] + field.Vpx.e[i][j][k]) * (field.Sxx.e[i][j][k] + field.Sxx.e[im][j][k]);
				VxSxx += (field.Vx.e[ip][j][k] + field.Vpx.e[ip][j][k]) * (field.Sxx.e[i][j][k] + field.Sxx.e[ip][j][k]);
				VxSxx *= -0.125;;
				// VySxy
				VySxy = field.Sxy.e[i][j][k] * (field.Vy.e[im][j][k] + field.Vy.e[i][j][k] + field.Vpy.e[im][j][k] + field.Vpy.e[i][j][k]);
				VySxy += field.Sxy.e[i][jp][k] * (field.Vy.e[im][jp][k] + field.Vy.e[i][jp][k] + field.Vpy.e[im][jp][k] + field.Vpy.e[i][jp][k]);
				VySxy += field.Sxy.e[ip][j][k] * (field.Vy.e[i][j][k] + field.Vy.e[ip][j][k] + field.Vpy.e[i][j][k] + field.Vpy.e[ip][j][k]);
				VySxy += field.Sxy.e[ip][jp][k] * (field.Vy.e[i][jp][k] + field.Vy.e[ip][jp][k] + field.Vpy.e[i][jp][k] + field.Vpy.e[ip][jp][k]);
				VySxy *= -0.0625;
				// VzSxz
				VzSxz = field.Sxz.e[i][j][k] * (field.Vz.e[im][j][k] + field.Vz.e[i][j][k] + field.Vpz.e[im][j][k] + field.Vpz.e[i][j][k]);
				VzSxz += field.Sxz.e[i][j][kp] * (field.Vz.e[im][j][kp] + field.Vz.e[i][j][kp] + field.Vpz.e[im][j][kp] + field.Vpz.e[i][j][kp]);
				VzSxz += field.Sxz.e[ip][j][k] * (field.Vz.e[i][j][k] + field.Vz.e[ip][j][k] + field.Vpz.e[i][j][k] + field.Vpz.e[ip][j][k]);
				VzSxz += field.Sxz.e[ip][j][kp] * (field.Vz.e[i][j][kp] + field.Vz.e[ip][j][kp] + field.Vpz.e[i][j][kp] + field.Vpz.e[ip][j][kp]);
				VzSxz *= -0.0625;
				// Jx
				jx += (VxSxx + VySxy + VzSxz);

				// y direction
				// VxSxy
				VxSxy = field.Sxy.e[i][j][k] * (field.Vx.e[i][jm][k] + field.Vx.e[i][j][k] + field.Vpx.e[i][jm][k] + field.Vpx.e[i][j][k]);
				VxSxy += field.Sxy.e[ip][j][k] * (field.Vx.e[ip][jm][k] + field.Vx.e[ip][j][k] + field.Vpx.e[ip][jm][k] + field.Vpx.e[ip][j][k]);
				VxSxy += field.Sxy.e[i][jp][k] * (field.Vx.e[i][j][k] + field.Vx.e[i][jp][k] + field.Vpx.e[i][j][k] + field.Vpx.e[i][jp][k]);
				VxSxy += field.Sxy.e[ip][jp][k] * (field.Vx.e[ip][j][k] + field.Vx.e[ip][jp][k] + field.Vpx.e[ip][j][k] + field.Vpx.e[ip][jp][k]);
				VxSxy *= -0.0625;
				// VySyy
				VySyy = (field.Vy.e[i][j][k] + field.Vpy.e[i][j][k]) * (field.Syy.e[i][j][k] + field.Syy.e[i][jm][k]);
				VySyy += (field.Vy.e[i][jp][k] + field.Vpy.e[i][jp][k]) * (field.Syy.e[i][j][k] + field.Syy.e[i][jp][k]);
				VySyy *= -0.125;
				// VzSyz
				VzSyz = field.Syz.e[i][j][k] * (field.Vz.e[i][jm][k] + field.Vz.e[i][j][k] + field.Vpz.e[i][jm][k] + field.Vpz.e[i][j][k]);
				VzSyz += field.Syz.e[i][j][kp] * (field.Vz.e[i][jm][kp] + field.Vz.e[i][j][kp] + field.Vpz.e[i][jm][kp] + field.Vpz.e[i][j][kp]);
				VzSyz += field.Syz.e[i][jp][k] * (field.Vz.e[i][j][k] + field.Vz.e[i][jp][k] + field.Vpz.e[i][j][k] + field.Vpz.e[i][jp][k]);
				VzSyz += field.Syz.e[i][jp][kp] * (field.Vz.e[i][j][kp] + field.Vz.e[i][jp][kp] + field.Vpz.e[i][j][kp] + field.Vpz.e[i][jp][kp]);
				VzSyz *= -0.0625;
				// jy
				jy += (VxSxy + VySyy + VzSyz);

				// z direction
				// VxSxz
				VxSxz = field.Sxz.e[i][j][k] * (field.Vx.e[i][j][km] + field.Vx.e[i][j][k] + field.Vpx.e[i][j][km] + field.Vpx.e[i][j][k]);
				VxSxz += field.Sxz.e[ip][j][k] * (field.Vx.e[ip][j][km] + field.Vx.e[ip][j][k] + field.Vpx.e[ip][j][km] + field.Vpx.e[ip][j][k]);
				VxSxz += field.Sxz.e[i][j][kp] * (field.Vx.e[i][j][k] + field.Vx.e[i][j][kp] + field.Vpx.e[i][j][k] + field.Vpx.e[i][j][kp]);
				VxSxz += field.Sxz.e[ip][j][kp] * (field.Vx.e[ip][j][k] + field.Vx.e[ip][j][kp] + field.Vpx.e[ip][j][k] + field.Vpx.e[ip][j][kp]);
				VxSxz *= -0.0625;
				// VySyz
				VySyz = field.Syz.e[i][j][k] * (field.Vy.e[i][j][km] + field.Vy.e[i][j][k] + field.Vpy.e[i][j][km] + field.Vpy.e[i][j][k]);
				VySyz += field.Syz.e[i][jp][k] * (field.Vy.e[i][jp][km] + field.Vy.e[i][jp][k] + field.Vpy.e[i][jp][km] + field.Vpy.e[i][jp][k]);
				VySyz += field.Syz.e[i][j][kp] * (field.Vy.e[i][j][k] + field.Vy.e[i][j][kp] + field.Vpy.e[i][j][k] + field.Vpy.e[i][j][kp]);
				VySyz += field.Syz.e[i][jp][kp] * (field.Vy.e[i][jp][k] + field.Vy.e[i][jp][kp] + field.Vpy.e[i][jp][k] + field.Vpy.e[i][jp][kp]);
				VySyz *= -0.0625;
				// VzSzz
				VzSzz = (field.Vz.e[i][j][k] + field.Vpz.e[i][j][k]) * (field.Szz.e[i][j][k] + field.Szz.e[i][j][km]);
				VzSzz += (field.Vz.e[i][j][kp] + field.Vpz.e[i][j][kp]) * (field.Szz.e[i][j][k] + field.Szz.e[i][j][kp]);
				VzSzz *= -0.125;
				// jz
				jz += (VxSxz + VySyz + VzSzz);
			}

	i = sim->energy.counter - 1;
	jx += sim->evx;
	jy += sim->evy;
	jz += sim->evz;

	i = sim->flux.counter;
	sim->flux.x[i] = jx;
	sim->flux.y[i] = jy;
	sim->flux.z[i] = jz;
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

	// normal exit
	return 0;
}
