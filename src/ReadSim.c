// This file reads the simulation setting
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include "Datastruct.h"
#include"functions.h"
#include<time.h>

int ReadSim(Sim * sim)
{
	FILE *fp;
	unsigned long long RandSeed;

	fp = fopen("Simulation.in","a+");
	// Number of times steps.
	fscanf(fp,"%d", &sim->Neq);
	fscanf(fp,"%d", &sim->Npr);
	fscanf(fp,"%d", &sim->Ns);			// sampling interval.
	// time step
	fscanf(fp,"%lf", &sim->dt);
	// grid spacing
	fscanf(fp,"%lf", &sim->dx);
	// elastic constants
	fscanf(fp,"%lf", &sim->lambda);
	fscanf(fp,"%lf", &sim->mu);
	fscanf(fp,"%lf", &sim->rho);
	// three compliance constant for strain energy calculation
	fscanf(fp,"%lf", &sim->s11);
	fscanf(fp,"%lf", &sim->s12);
	fscanf(fp,"%lf", &sim->s44);
	fclose(fp);

	// calculate the velocities and the maximum allowable time step for FDTD stability
	sim->lam2mu = sim->lambda + 2*sim->mu;
	sim->vp = sqrt(sim->lam2mu / sim->rho);
	sim->vs = sqrt(sim->mu / sim->rho);
	sim->dV = sim->dx * sim->dx * sim->dx;
	sim->mass = sim->dV * sim->rho;
	sim->mdt = sim->dx / sqrt(3) / sim->vp;

	// Hard Code Output interval
	sim->No = 100000;

	// print out to the screen
	printf("===========Overview of the Simulation=============\n");
	printf("Equilibration: %d steps\n", sim->Neq);
	printf("Production: %d steps\n", sim->Npr);
	printf("Sampling interval: %d steps\n", sim->Ns);
	printf("Time step: %e (s)\n", sim->dt);
	printf("Grid spacing: %e (m)\n", sim->dx);
	printf("Lambda: %e (Pa)\n", sim->lambda);
	printf("Mu: %e (Pa)\n", sim->mu);
	printf("Mass density: %e (kg/m^3)\n", sim->rho);
	printf("Longitudinal Velocity: %e (m/s)\n", sim->vp);
	printf("Transverse Velocity: %e (m/s)\n", sim->vs);
	printf("CFL time step: %e (s)\n", sim->mdt);

	// output the simulation setting
	// Initialize the RNG
	// Initialize RNG seed
	srand((unsigned)time(NULL));
	RandSeed = (unsigned long long)((rand() * 9973 + 9949) * 7211 + 6067);
	InitRand(RandSeed);
	printf("Rand Seed:%llu\n", RandSeed);

	fp = fopen("Simulation.log","a+");
	fprintf(fp,"=======Simulation Setting==========\n");
	fprintf(fp,"Rand Seed: %llu\n", RandSeed);
	fprintf(fp,"Num. of Equilibration Steps: %d\n", sim->Neq);
	fprintf(fp,"Num. of Production Steps: %d\n", sim->Npr);
	fprintf(fp,"Sampling interval: %d steps\n", sim->Ns);
	fprintf(fp,"Output interval: %d samplings\n", sim->No);
	fprintf(fp,"Grid spacing:%e (m)\n", sim->dx);
	fprintf(fp,"Time step:%e (s)\n", sim->dt);
	fprintf(fp,"Lambda: %e (Pa)\n", sim->lambda);
	fprintf(fp,"Mu: %e (Pa)\n", sim->mu);
	fprintf(fp,"Vp: %e (m/s) \n", sim->vp);
	fprintf(fp,"Vs: %e (m/s)\n", sim->vs);
	fprintf(fp, "s11: %e\n", sim->s11);
	fprintf(fp, "s12: %e\n", sim->s12);
	fprintf(fp, "s44: %e\n", sim->s44);
	fclose(fp);
	// normal exit
	return 0;
}
