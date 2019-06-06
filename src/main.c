// The controlling over the flow
#include"Datastruct.h"
#include"functions.h"

// the main controlling flow function.
int main()
{
	Sim sim;							// global simulation setting
	Field field;						// simulation field
	Coeff coeff;						// interpolated coefficient field
	Mat rho;							// simulation domain

	int ti;

	// read in the global simulation setting.
	ReadSim(&sim);
	// Need to allocate memory for the simulation field and coefficient matrices
	// read the density matrix (the simulaiton domain)
	ReadDensityMat(&rho);
	// Initialize simulation

	InitField(&sim, rho, &field, &coeff);

	// Equilibration Stage.
	for (ti = 0; ti < sim.Neq; ti++)
	{
		UpdateStress(sim, field, coeff);
		GlobalEnergy(sim, rho, field);
		UpdateVelocity(sim, field, coeff);
	}

	// Production Stage
	for (ti = 0; ti < sim.Npr; ti++)
	{
		UpdateStress(sim, field, coeff);
		GlobalEnergy(sim, rho, field);
		HeatCurrent(sim, rho, field);
		UpdateVelocity(sim, field, coeff);
	}
	// normal exit
	return 0;
}
