/*
 * singlePendulumSim.cpp
 *
 *  Created on: Feb 24, 2023
 *      Author: norms
 */

#include "../headers/singlePendulumSim.h"
#include "../headers/constraint_bodies.h"
#include "../headers/get_state.h"
#include "../headers/utils.h"
#include "../headers/renderer.h"
#include "../headers/matrix_stuff.h"
#include "../headers/rk4.h"
#include "../headers/rigid_bodies.h"
#include "../headers/springs.h"
#include "../headers/toFile.h"
#include "../headers/total_energy.h"

// Simulates a single pendulum of specified customization and outputs data to file
void singlePendulumSim(int numMasses, Circular_Rigid_Body *masses, int numBar1s, Rigid_Bar_1 *bar1s,  int numBar2s, Rigid_Bar_2 *bar2s, double dt, int totalTime)
{
	// Initialize the springs: There arent actually any
	Spring_2 spring2s[1];

	// Create the State_Getter class for current state
	State_Getter current_state;

	// Number of masses and constraints
	current_state.num_bodies = numMasses;
	current_state.num_constraints = numBar1s + numBar2s;

	// Total Simulation Timesteps
	int totalSteps = totalTime / dt;

	// Data collection
	// Data array for storage of position data: Store the pivot position, and mass positions

	// Collect data every ITERATIONS step
	int ITERATIONS = 1000;
	int dataLen = ((totalSteps + ITERATIONS - 1) / ITERATIONS) + 1;

	// Data for pivot point
	Vector pivotData = {bar1s[0].pivot.x, bar1s[0].pivot.y};

	// Data for Masses
	Vector **massData = (Vector **)malloc(dataLen * sizeof(Vector *));

	// Data for angles
	int numAngles = numMasses + 1;

	double **angleData = (double **)malloc(dataLen * sizeof(double *));

	// Allocate space for the masses
	for (int i = 0; i < dataLen; ++i)
	{
		massData[i] = (Vector *)malloc(numMasses * sizeof(Vector));
		angleData[i] = (double *)malloc(numAngles * sizeof(double));
	}

	// Loop
	// Loop count
	int loopCount = 0;
	// Data bucket count
	int dc = 0;
	while (loopCount < totalSteps)
	{
		// Store data
		if (loopCount % ITERATIONS == 0)
		{
			// Masses
			for (int i = 0; i < numMasses; ++i)
			{
				massData[dc][i] = {masses[i].pos.x, masses[i].pos.y};
			}

			// Angles
			// Angle between pivot and first mass
			Vector vPivot = {masses[0].pos.x - bar1s[0].pivot.x, masses[0].pos.y - bar1s[0].pivot.y};
			Vector vert = {0, -1};
			double angle1 = angle_between_vectors(vPivot, vert);
			angleData[dc][0] = angle1;

			// Find the rest of the angles
			Vector v;
			double angle;
			if (numAngles > 1)
			{
				for (int i = 1; i < numMasses; ++i)
				{
					v = {masses[i].pos.x - masses[i - 1].pos.x, masses[i].pos.y - masses[i - 1].pos.y};
					angle = angle_between_vectors(v, vert);
					angleData[dc][i] = angle;
				}
			}

//			std::cout << "Angles" << std::endl;
//			std::cout << "Angle 1: " << angleData[dc][0] * 180 / M_PI << std::endl;
//			std::cout << "Angle 2: " << angleData[dc][1] * 180 / M_PI << std::endl;
//			std::cout << std::endl;
//			std::cout << "Constraints" << std::endl;
//			std::cout << "c1: " << bar1s[0].constraint(masses) << std::endl;
//			std::cout << "c2: " << bar2s[0].constraint(masses) << std::endl;
//			std::cout << "Total energy: " << total_energy(masses, numMasses, loopCount) << std::endl << std::endl;

			// Increment data bucket count
			dc++;

		}

		// Numerical integrator
		// rk4 method
		rk4(&current_state, masses, bar1s, numBar1s, bar2s, numBar2s, spring2s, dt);

		// Increment the loop count
		loopCount++;
	}

	// Record last bucket of data for masses
	for (int i = 0; i < numMasses; ++i)
	{
		massData[dataLen - 1][i] = {masses[i].pos.x, masses[i].pos.y};
	}

	// Record final angles
	// Angle between pivot and first mass
	Vector vPivot = {masses[0].pos.x - bar1s[0].pivot.x, masses[0].pos.y - bar1s[0].pivot.y};
	Vector vert = {0, -1};
	double angle1 = angle_between_vectors(vPivot, vert);
	angleData[dataLen - 1][0] = angle1;

	// Find the rest of the angles
	Vector v;
	double angle;
	if (numAngles > 1)
	{
		for (int i = 1; i < numMasses; ++i)
		{
			v = {masses[i].pos.x - masses[i - 1].pos.x, masses[i].pos.y - masses[i - 1].pos.y};
			angle = angle_between_vectors(v, vert);
			angleData[dataLen - 1][i] = angle;
		}
	}

	// Write data to file
	pendulumToFile(pivotData, massData, dataLen, numMasses);

	// Success
	std::cout << "Done" << std::endl;
}
