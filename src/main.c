// The controlling over the flow
#include"Datastruct.h"
#include"functions.h"
#include<stdio.h>
#include<time.h>

// the main controlling flow function.
int main()
{
	time_t start, end;
	Sim sim;							// global simulation setting
	Field field;						// simulation field
	Coeff coeff;						// interpolated coefficient field
	Mat rho;							// simulation domain
	int ti;
	FILE *fp;

	time(&start);

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
		UpdateStress(sim, rho, field, coeff);
		UpdateVelocity(sim, field, coeff);

		if (ti % sim.Ns == 0)
			SampleEnergy(&sim, rho, field);
	}

	// Production Stage
	for (ti = 0; ti < sim.Npr; ti++)
	{
		UpdateStress(sim, rho, field, coeff);
		UpdateVelocity(sim, field, coeff);

		if (ti % sim.Ns == 0)
		{
			SampleEnergy(&sim, rho, field);
			SampleHeatCurrent(&sim, rho, field);
		}
	}
	time(&end);

	fp = fopen("Simulation.log","a+");
	fprintf(fp,"Elapsed time: %.2f (s)\n", difftime(end, start));
	fprintf(fp,"EiSumAve: %e\n", sim.EiSumAve * sim.time / (sim.Ns * sim.dt));
	fclose(fp);
	// normal exit
	return 0;
}
