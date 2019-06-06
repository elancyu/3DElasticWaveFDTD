// This file contains the declaration of all functions.
#include "Datastruct.h"

/*==================================Energy.c===========================================*/
// Calculate the global energy
int GlobalEnergy(Sim sim, Mat rho,Field field);

// Calculate the heat current
int HeatCurrent(Sim sim, Mat rho, Field field);

/*==================================Init.c===========================================*/
// initialize simulation field
int InitField(Sim *sim, Mat rho, Field *field, Coeff *coeff);

// Initialize velocity
int InitVelocity(Sim *sim, Mat rho, Field *field, Coeff *coeff);

// Interpolate coefficients including buoyancy and elastic constants.
int InterpolateCoeffs(Mat rho, Coeff *coeff);

// Multiply the coefficient matrix with the actual physical parameters
int PhysicalCoeff(Sim sim, Coeff *coeff);

/*==================================Matrix.c===========================================*/
// set the 3D matrix size
int SetMatSize(Mat *M, int Nx, int Ny, int Nz);

// Allocate memory and initialze for the 3D matrix.
int MallocMat(Mat *M, double InitVal);

/*==================================Rand.c===========================================*/
// initialize random number generator with assigned random seed.
void InitRand(unsigned long long seed);

// Generate a random number uniformly distributed within (0,1)
double RandR(void);

// Gernate a random number following normal distirubtion using Box-Muller algorithm, E is expectation and V is standard variance.
double RandG(double E, double V);

/*==================================ReadGeo.c===========================================*/
// read in the geometry, 1 for solid and 0 for air.
int ReadDensityMat(Mat *rho);

/*==================================ReadGeo.c===========================================*/
int ReadDensityMat(Mat *rho);

/*==================================ReadSim.c===========================================*/
// read in the simulation settings and do the preliminary initialization
int ReadSim(Sim * sim);

/*==================================UpdateStress.c===========================================*/
// update the stress field.
int UpdateStress(Sim sim, Field field, Coeff coeff);

/*==================================UpdateVelocity.c===========================================*/
// update the velocity field
int UpdateVelocity(Sim sim, Field field, Coeff coeff);
