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
	Mat Vx;
	Mat Vy;
	Mat Vz;
	Mat Sxx;
	Mat Syy;
	Mat Szz;
	Mat Syz;
	Mat Sxz;
	Mat Sxy;
}Field;

// the coefficient matrices
typedef struct Coeff{
	Mat buox;
	Mat buoy;
	Mat buoz;
	Mat Muyz;
	Mat Muxz;
	Mat Muxy;
}Coeff;

#endif
