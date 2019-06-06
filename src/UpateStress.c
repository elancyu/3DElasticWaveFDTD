// This file updates the stress, assuming periodic boundary conditions in X and Y directions.

#include "Datastruct.h"
#include<math.h>

int UpdateStress(Sim sim, Field field, Coeff coeff)
{
	double c11, c12, c44;
	double dVxdx, dVydy, dVzdz;
	double dtdx;				// dt/dx
	int i, j, k;
	int im, ip, jm, jp, km, kp;			// the neighboring indices.
	int Nx, Ny, Nz;

	Nx = field.Nx;
	Ny = field.Ny;
	Nz = field.Nz;
	dtdx = sim.dt / sim.dx;

	c11 = 0.0;
	c12 = 0.0;
	c44 = 0.0;

	// update the simulation time
	sim.time += 0.5 * sim.dt;

	// Maybe better way to deal with the updating
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 1; k < Nz - 1; k++)
			{
				// first calculate the divergence
				ip = (i + 1) % Nx;
				jp = (j + 1) % Ny;
				kp = k + 1;
				dVxdx = field.Vx.e[ip][j][k] - field.Vx.e[i][j][k];
				dVydy = field.Vy.e[i][jp][k] - field.Vy.e[i][j][k];
				dVzdz = field.Vz.e[i][j][kp] - field.Vz.e[i][j][k];
				// Update the normal stress
				field.Sxx.e[i][j][k] += dtdx * (coeff.Lam2mu.e[i][j][k] * dVxdx + coeff.Lam.e[i][j][k] * dVydy + coeff.Lam.e[i][j][k] * dVzdz);
				field.Syy.e[i][j][k] += dtdx * (coeff.Lam.e[i][j][k] * dVxdx + coeff.Lam2mu.e[i][j][k] * dVydy + coeff.Lam.e[i][j][k] * dVzdz);
				field.Szz.e[i][j][k] += dtdx * (coeff.Lam.e[i][j][k] * dVxdx + coeff.Lam.e[i][j][k] * dVydy + coeff.Lam2mu.e[i][j][k] * dVzdz);

				// Update the shear stress
				im = (i - 1 + Nx) % Nx;
				jm = (j - 1 + Ny) % Ny;
				field.Sxy.e[i][j][k] += coeff.Muxy.e[i][j][k] * dtdx * (field.Vy.e[i][j][k] - field.Vy.e[im][j][k] + field.Vx.e[i][j][k] - field.Vx.e[i][jm][k]);
			}

	// Update the other two shear stress
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 1; k < Nz; k++)
			{
				im = (i - 1 + Nx) % Nx;
				jm = (j - 1 + Ny) % Ny;
				km = k - 1;
				field.Sxz.e[i][j][k] += coeff.Muxz.e[i][j][k] * dtdx * (field.Vz.e[i][j][k] - field.Vz.e[im][j][k] + field.Vx.e[i][j][k] - field.Vx.e[i][j][km]);
				field.Syz.e[i][j][k] += coeff.Muyz.e[i][j][k] * dtdx * (field.Vz.e[i][j][k] - field.Vz.e[i][jm][k] + field.Vy.e[i][j][k] - field.Vy.e[i][j][km]);
			}
	// normal exit
	return 0;
}
