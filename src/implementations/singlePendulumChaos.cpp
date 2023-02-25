/*
 * singlePendulumChaos.cpp
 *
 *  Created on: Feb 24, 2023
 *      Author: norms
 */


#include "../headers/singlePendulumChaos.h"
#include "../headers/circleGenerator.h"

// Simulates a single double pendulum with a circle of perturbations in the angle space
void singleDoublePendulumChaos(double *angles, double *barLens, double *masses, double dt, int totalTime, int numPoints, double radius)
{
	// Write the unit circle of perturbations to file
	writeCircleToFile(numPoints, radius);

	// Initialize the pendulum

	// Initialize perturbations of the pendulum



}

