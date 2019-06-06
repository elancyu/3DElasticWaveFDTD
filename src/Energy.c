// This file contains the energy, momentum and flux calculations.
#include"functions.h"
#include"Datastruct.h"
#include<stdio.h>

// Strictly calculate the kinetic energy and strain energy.
int GlobalEnergy(Sim sim, Mat rho,Field field)
{
	FILE *fp;
	int i, j, k;
	int im, jm, km;
	double wyz, wxz, wxy;
	double wx, wy, wz;
	double Ek, Es;
	int Nx, Ny, Nz;
	// remember to remove s11, s12, s44;
	Es = 0.0;
	Ek = 0.0;

	// retrieve the matrix size
	Nx = rho.Nx;
	Ny = rho.Ny;
	Nz = rho.Nz;

	// directly add up all Sii^2, and the cross term, SiiSjj.
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
			{
				// minus direction indices
				im = (i - 1 + Nx) % Nx;
				jm = (j - 1 + Ny) % Ny;
				km = (k - 1 + Nz) % Nz;
				// Sii^2 && SiiSjj
				Es += 0.5 * sim.s11 * rho.e[i][j][k] * field.Sxx.e[i][j][k] * field.Sxx.e[i][j][k];
				Es += 0.5 * sim.s11 * rho.e[i][j][k] * field.Syy.e[i][j][k] * field.Syy.e[i][j][k];
				Es += 0.5 * sim.s11 * rho.e[i][j][k] * field.Szz.e[i][j][k] * field.Szz.e[i][j][k];
				Es += sim.s12 * rho.e[i][j][k] * field.Sxx.e[i][j][k] * field.Syy.e[i][j][k];
				Es += sim.s12 * rho.e[i][j][k] * field.Sxx.e[i][j][k] * field.Szz.e[i][j][k];
				Es += sim.s12 * rho.e[i][j][k] * field.Syy.e[i][j][k] * field.Szz.e[i][j][k];
				// weighting factor for Sij^2
				wyz = 0.25 * (rho.e[i][j][k] + rho.e[i][jm][k] + rho.e[i][j][km] + rho.e[i][jm][km]);
				wxz = 0.25 * (rho.e[i][j][k] + rho.e[im][j][k] + rho.e[i][j][km] + rho.e[im][j][km]);
				wxy = 0.25 * (rho.e[i][j][k] + rho.e[im][j][k] + rho.e[i][jm][k] + rho.e[im][jm][k]);
				// Sij^2
				Es += sim.s44 * wyz * field.Syz.e[i][j][k] * field.Syz.e[i][j][k];
				Es += sim.s44 * wxz * field.Sxz.e[i][j][k] * field.Sxz.e[i][j][k];
				Es += sim.s44 * wxy * field.Sxy.e[i][j][k] * field.Sxy.e[i][j][k];

				// Kinetic Energy
				wx = 0.5 * (rho.e[i][j][k] + rho.e[im][j][k]);
				wy = 0.5 * (rho.e[i][j][k] + rho.e[i][jm][k]);
				wz = 0.5 * (rho.e[i][j][k] + rho.e[i][j][km]);
				Ek += 0.5 * wx * sim.rho * field.Vx.e[i][j][k] * field.Vx.e[i][j][k];
				Ek += 0.5 * wy * sim.rho * field.Vy.e[i][j][k] * field.Vy.e[i][j][k];
				Ek += 0.5 *wz * sim.rho * field.Vz.e[i][j][k] * field.Vz.e[i][j][k];
			}

	// times the volume
	Es *= sim.dV;
	Ek *= sim.dV;

	// output the energy
	fp = fopen("Energy.log","a+");
	fprintf(fp,"%e  %e  %e\n", Ek, Es, Ek + Es);
	fclose(fp);
	// normal exit
	return 0;
}

// Sample the heat current spatial integral, first interpolate to velocity and then integrate over the domain
int HeatCurrent(Sim sim, Mat rho, Field field)
{
	FILE *fp;
	int i, j, k;
	int im, jm, km, ip, jp, kp;			// interpolating indices;
	int Nx, Ny, Nz;
	double wx, wy, wz;					// weighting factor for integration.
	double VxSxx, VySxy, VzSxz; 		// components of Jx
	double VxSxy, VySyy, VzSyz;			// components of Jy
	double VxSxz, VySyz, VzSzz;			// components of Jz
	double jx, jy, jz;

	Nx = rho.Nx;
	Ny = rho.Ny;
	Nz = rho.Nz;

	jx = 0.0;
	jy = 0.0;
	jz = 0.0;

	// quantity related to Vx and Vy
	for (i = 0; i < Nx; i++)
		for (j = 0; j < Ny; j++)
			for (k = 0; k < Nz; k++)
			{
				// interpolating indices
				im = (i - 1 + Nx) % Nx;
				jm = (j - 1 + Ny) % Ny;
				km = (k - 1 + Nz) % Nz;
				ip = (i + 1) % Nx;
				jp = (j + 1) % Ny;
				kp = (k + 1) % Nz;
				// Jx: VxSxx, interpolation
				VxSxx = field.Vx.e[i][j][k] * (field.Sxx.e[i][j][k] + field.Sxx.e[im][j][k]);
				// Jx: VySxy
				VySxy = field.Vy.e[i][j][k] * (field.Sxy.e[i][j][k] + field.Sxy.e[i][jp][k]);
				// Jx: VzSxz, interpolation
				VzSxz = field.Vz.e[i][j][k] * (field.Sxz.e[i][j][k] + field.Sxz.e[i][j][kp]);
				// Jy: VxSxy, interpolation
				VxSxy = field.Vx.e[i][j][k] * (field.Sxy.e[i][j][k] + field.Sxy.e[ip][j][k]);
				// Jy: VySyy, interpolation
				VySyy = field.Vy.e[i][j][k] * (field.Syy.e[i][j][k] + field.Syy.e[i][jm][k]);
				// jy: VzSyz
				VzSyz = field.Vy.e[i][j][k] * (field.Syz.e[i][j][k] + field.Syz.e[i][j][kp]);
				// Jz: VxSxz, interpolation
				VxSxz = field.Vx.e[i][j][k] * (field.Sxz.e[i][j][k] + field.Sxz.e[ip][j][k]);
				// Jz: VySyz, interpolation
				VySyz = field.Vy.e[i][j][k] * (field.Syz.e[i][j][k] + field.Syz.e[i][jp][k]);
				// Jz: VzSzz
				VzSzz = field.Vz.e[i][j][k] * (field.Szz.e[i][j][k] + field.Szz.e[i][j][km]);

				// weighting factor
				wx = 0.5 * (rho.e[i][j][k] + rho.e[im][j][k]);
				wy = 0.5 * (rho.e[i][j][k] + rho.e[i][jm][k]);
				wz = 0.5 * (rho.e[i][j][k] + rho.e[i][j][km]);
				// Trapezoidal Integration
				// X direction: related to Vx
				jx += wx * VxSxx;
				jy += wx * VxSxy;
				jz += wx * VxSxz;
				// Y direction: related to Vy
				jx += wy * VySxy;
				jy += wy * VySyy;
				jz += wy * VySyz;
				// Z direction: related to Vz
				jx += wz * VzSxz;
				jy += wz * VzSyz;
				jz += wz * VzSzz;
			}
	// output the heat current integral to file
	fp = fopen("HC.out","a+");
	fprintf(fp,"%e  %e  %e\n", jx, jy, jz);
	fclose(fp);
	// normal exit
	return 0;
}
