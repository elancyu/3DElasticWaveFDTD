// This file defines the simulation setting
#ifndef DATASTRUCT_H
#define DATASTRUCT_H

// the data sturcture for the energy
typedef struct Energy{
	int len;		// indicate the length of the array, do not change it.
	int counter;	// count the valid data number.
	double *ke;
	double *se;
	double *te;
}Energy;

// the data structure for the heat current
typedef struct Flux{
	int len;			// indicate the length of the array, do not change it.
	int counter;		// count the valid data number.
	double *x;
	double *y;
	double *z;
}Flux;

// definition of the 3D matrix
typedef struct Mat{
	int Nx;				// X size
	int Ny;				// Y size
	int Nz;				// Z size
	double ***e;		// element matrix
}Mat;

// definition of simulation settings
typedef struct Sim{
	int Neq;						// number of time steps in equilibration stage.
	int Npr;						// number of time steps in production stage.
	int Ns;							// Sampling interval to reduce I/O time.
	double dx;						// grid spacing.
	double dt;						// time step.
	double lambda;					// lame constant lambda
	double lam2mu;					// lambda + 2mu in the cell center.
	double mu;						// lame constant mu
	double rho;						// mass density.
	double vp;						// longitudinal wave speed.
	double vs;						// transverse wave speed.
	double mdt;						// maximum allowable time step.
	double dV;						// volume for each cell.
	double mass;					// mass for each cell.
	int Natom;						// number of nozero mass density for the cells.
	double time;					// simulaiton time recorder
	double s11, s12, s44;			// simply convenient for the strain energy calculation.
	struct Energy energy;			// energy records
	struct Flux flux;				// flux records.
}Sim;

// the simulation field.
typedef struct Field{
	int Nx, Ny, Nz;					// size of the field.
	Mat Vx;							// the velocity at x direction, used for computation
	Mat Vy;							// the velocity at y direciton, used for computation
	Mat Vz;							// the velocity at z direction, used for computation
	Mat Vpx;						// previous step velocity at x direction.
	Mat Vpy;						// previous step velocity at y direciton.
	Mat Vpz;						// previous step velocity at z direction.
	Mat Sxx;						// normal stress 1
	Mat Syy;						// normal stress 2
	Mat Szz;						// normal stress 3
	Mat Syz;						// shear stress 1
	Mat Sxz;						// shear stress 2
	Mat Sxy;						// sheat stress 3
}Field;

// the coefficient matrices
typedef struct Coeff{
	Mat buox;						// buoyancy at x
	Mat buoy;						// buoyancy at y
	Mat buoz;						// buoyancy at z
	Mat Muyz;						// interpolated elastic constant at yz.
	Mat Muxz;						// interpolated elastic constant at xz.
	Mat Muxy;						// interpolated elastic constant at xy.
}Coeff;

#endif
