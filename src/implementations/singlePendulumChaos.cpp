/*
 * singlePendulumChaos.cpp
 *
 *  Created on: Feb 24, 2023
 *      Author: norms
 */


#include "../headers/singlePendulumChaos.h"
#include "../headers/circleGenerator.h"
#include "../headers/Pendulum.h"

// Simulates a single double pendulum with a circle of perturbations in the angle space
void singleDoublePendulumChaos(double *angles, double *barLens, double *massAmounts, double dt, int totalTime, int numPoints, double radius)
{
	// Write the unit circle of perturbations to file
	writeCircleToFile(numPoints, radius);

	// Initialize the double pendulum
	Pendulum *p = initializeDoublePendulum(angles, barLens, massAmounts);

	// Initialize perturbations of the pendulum


}

// Initializes perturbations of a double pendulum using offsets from file
void initializePerturbations(Pendulum *p, int numPoints)
{
	// Allocate space for the perturbations
	p->perturbed = (Pendulum *)malloc(numPoints * sizeof(Pendulum));

	// Copy pendulum data into each perturbation
}

