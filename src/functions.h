// This file contains the declaration of all functions.
#include "Datastruct.h"

/*==================================Energy.c===========================================*/
// Calculate the global energy
int SampleEnergy(Sim *sim, Mat rho,Field field);

// Calculate the heat current
int SampleHeatCurrent(Sim *sim, Mat rho, Field field);

// write out the energy and then reset the data
int WriteEnergy(Sim *sim);

// write out the heat current and reset the data
int WriteHeatCurrent(Sim *sim);

/*==================================Init.c===========================================*/
// initialize simulation field
int InitField(Sim *sim, Mat rho, Field *field, Coeff *coeff);

// Initialize velocity
int InitVelocity(Sim *sim, Mat rho, Field *field, Coeff *coeff);

// Interpolate coefficients including buoyancy and elastic constants, for free boundary condition implementation
int InterpolateFreeCoeff(Mat rho, Coeff *coeff);

// Interpolate coefficients including buoyancy and elastic constants, for fixed boundary condition implementation
int InterpolateFixedCoeff(Mat rho, Coeff *coeff);

// Multiply the coefficient matrix with the actual physical parameters
int PhysicalCoeff(Sim sim, Coeff *coeff);

// Initialize the sampling recorders
int InitSampling(Sim *sim, int len);

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
int UpdateStress(Sim sim, Mat rho, Field field, Coeff coeff);

/*==================================UpdateVelocity.c===========================================*/
// update the velocity field
int UpdateVelocity(Sim sim, Field field, Coeff coeff);

/*==================================NormalTest.c===========================================*/
// Calculate the Jarque-Bera value for the distribution.
double JBTest(double *dist, int len);

// Generate a random number array in normal distribution so that the JB value is lower than some threshold.
int GenRandGArray(double *arr, int len, int Natom);

// Initialize the velocity in 3 directions.
int InitVelocityS(Sim *sim, Mat rho, Field *field, Coeff *coeff);
