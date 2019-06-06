// This file updates the velocity, assuming periodic boundary condition in X and Y directions.
// Note: time step and spacing grid is not included yet.

#include "Datastruct.h"
#include<math.h>

int UpdateVelocity(Sim sim, Field field, Coeff coeff)
{
	int i, j, k;
	int im, jm, km, ip, jp, kp;			// indices for differencing.
	int Nx, Ny, Nz;
	double dtdx;						// dt/dx

	Nx = field.Nx;
	Ny = field.Ny;
	Nz = field.Nz;
	dtdx = sim.dt/sim.dx;

	// update the simulation time
	sim.time += 0.5 * sim.dt;

	// update Vx and Vy
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 1; k < Nz - 1; k++)
			{
				// the neighbors for differencing.
				im = (i - 1 + Nx) % Nx;
				ip = (i + 1) % Nx;
				jm = (j - 1 + Ny) % Ny;
				jp = (j + 1) % Ny;
				kp = k + 1;
				// Vx
				field.Vx.e[i][j][k] += coeff.buox.e[i][j][k] * dtdx * (field.Sxx.e[i][j][k] - field.Sxx.e[im][j][k] + field.Sxy.e[i][jp][k] - field.Sxy.e[i][j][k] + field.Sxz.e[i][j][kp] - field.Sxz.e[i][j][k]);
				field.Vy.e[i][j][k] += coeff.buoy.e[i][j][k] * dtdx * (field.Sxy.e[ip][j][k] - field.Sxy.e[i][j][k] + field.Syy.e[i][j][k] - field.Syy.e[i][jm][k] + field.Syz.e[i][j][kp] - field.Syz.e[i][j][k]);
			}

	// Update Vz
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 1; k < Nz; k++)
			{
				ip = (i + 1) % Nx;
				jp = (j + 1) % Ny;
				km = k - 1;
				field.Vz.e[i][j][k] += coeff.buoz.e[i][j][k] * dtdx * (field.Sxz.e[ip][j][k] - field.Sxz.e[i][j][k] + field.Syz.e[i][jp][k] - field.Syz.e[i][j][k] + field.Szz.e[i][j][k] - field.Szz.e[i][j][km]);
			}
	// normal exit
	return 0;
}
