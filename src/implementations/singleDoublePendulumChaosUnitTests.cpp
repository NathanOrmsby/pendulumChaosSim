/*
 * singleDoublePendulumChaosUnitTests.cpp
 *
 *  Created on: Feb 24, 2023
 *      Author: norms
 */

#include <math.h>
#include <iostream>

#include "../headers/singleDoublePendulumChaosUnitTests.h"
#include "../headers/rigid_bodies.h"
#include "../headers/vectors.h"

// Unit test functions for single Double Pendulum Chaos simulation


// Unit test for initializeMassPositions Pendulum Function
Circular_Rigid_Body *initializeMassPositionsTest(double *angles, int numAngles, double *barLens, int numBars)
{
	std::cout << "\n\n\n" << "Unit test for initializing mass positions based on given angles and bar lengths\n" << std::endl;
	// Number of masses is equivalent to number of angles
	int numMasses = numAngles;

	// Allocate space
	Circular_Rigid_Body *masses = (Circular_Rigid_Body *)malloc(numMasses * sizeof(Circular_Rigid_Body));

	// Find cartesian coordinates for all masses
	// Running offset for shifted origin
	Vector offset = {0, 0};

	std::cout << "Calculating mass positions" << std::endl << std::endl;
	for (int i = 0; i < numMasses; ++i)
	{
		std::cout << "Mass: " << i << std::endl;

		masses[i].pos.x = barLens[i] * sin(angles[i]) + offset.x;
		masses[i].pos.y = barLens[i] * cos(angles[i]) + offset.y;

		std::cout << "Origin is: x: " << offset.x << " y: " << offset.y << std::endl;
		std::cout << "Angle is: " << angles[i] << std::endl;
		std::cout << "Bar length is: " << barLens[i] << std::endl;
		std::cout << "Mass position is: x: " << masses[i].pos.x << " y: " << masses[i].pos.y << std::endl;

		// Update the offset
		offset.x += masses[i].pos.x;
		offset.y += masses[i].pos.y;

		std::cout << "Offset is now: x: " << offset.x << " y: " << offset.y << std::endl;
	}
	std::cout << "\n\n-------------------------------------------------------------------------------------------------------------------------------------------------------\n\n\n";

	return masses;
}
