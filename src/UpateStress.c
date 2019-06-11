// This file updates the stress, assuming periodic boundary conditions in X and Y directions.
// Checked. The code is safe.

#include "Datastruct.h"
#include<math.h>

int UpdateStress(Sim sim, Mat rho, Field field, Coeff coeff)
{
	double c11, c12;
	double dVxdx, dVydy, dVzdz;
	double dtdx;				// dt/dx
	int i, j, k;
	int im, ip, jm, jp, km, kp;			// the neighboring indices.
	int Nx, Ny, Nz;

	Nx = field.Nx;
	Ny = field.Ny;
	Nz = field.Nz;
	dtdx = sim.dt / sim.dx;

	c11 = sim.lam2mu;
	c12 = sim.lambda;

	// update the simulation time
	sim.time += 0.5 * sim.dt;

	// Maybe better way to deal with the updating
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
			{
				// first calculate the divergence
				ip = (i + 1) % Nx;
				jp = (j + 1) % Ny;
				kp = (k + 1) % Nz;
				im = (i - 1 + Nx) % Nx;
				jm = (j - 1 + Ny) % Ny;
				km = (k - 1 + Nz) % Nz;

				// calculate divergence
				dVxdx = field.Vx.e[ip][j][k] - field.Vx.e[i][j][k];
				dVydy = field.Vy.e[i][jp][k] - field.Vy.e[i][j][k];
				dVzdz = field.Vz.e[i][j][kp] - field.Vz.e[i][j][k];

				// Update the normal stress
				field.Sxx.e[i][j][k] += dtdx * rho.e[i][j][k] * (c11 * dVxdx + c12 * dVydy + c12 * dVzdz);
				field.Syy.e[i][j][k] += dtdx * rho.e[i][j][k] * (c12 * dVxdx + c11 * dVydy + c12 * dVzdz);
				field.Szz.e[i][j][k] += dtdx * rho.e[i][j][k] * (c12 * dVxdx + c12 * dVydy + c11 * dVzdz);

				// Update the shear stress
				field.Sxy.e[i][j][k] += coeff.Muxy.e[i][j][k] * dtdx * (field.Vy.e[i][j][k] - field.Vy.e[im][j][k] + field.Vx.e[i][j][k] - field.Vx.e[i][jm][k]);
				field.Sxz.e[i][j][k] += coeff.Muxz.e[i][j][k] * dtdx * (field.Vz.e[i][j][k] - field.Vz.e[im][j][k] + field.Vx.e[i][j][k] - field.Vx.e[i][j][km]);
				field.Syz.e[i][j][k] += coeff.Muyz.e[i][j][k] * dtdx * (field.Vz.e[i][j][k] - field.Vz.e[i][jm][k] + field.Vy.e[i][j][k] - field.Vy.e[i][j][km]);
			}
	// normal exit
	return 0;
}
