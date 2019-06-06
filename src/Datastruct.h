// This file defines the simulation setting
#ifndef DATASTRUCT_H
#define DATASTRUCT_H
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
	double dx;						// grid spacing.
	double dt;						// time step.
	double lambda;					// lame constant lambda
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
	Mat Lam2mu;
	Mat Lam;
}Coeff;

#endif
