/*
 * singleDoublePendulumChaos.cpp
 *
 *  Created on: Apr 12, 2023
 *      Author: norms
 */

#include "../headers/Vectors.h"
#include "../headers/rk4.h"
#include "../headers/total_energy.h"
#include "../headers/DoublePendulum.h"
#include "../headers/singleDoublePendulumChaos.h"
#include "../headers/vectors.h"

#include <math.h>
#include <vector>
#include <tuple>
#include <iomanip>
#include <iostream>
#include <limits>
#include <algorithm>


// Timestep lyapunov function
void singleDoublePendulumTimestepLyapunov(Circular_Rigid_Body *masses, double *lyapunovRunningSums, int numPerturbations, double initialSeparation) {

	// Calculate angles of base pendulum
	std::vector<double> baseAngles = DoublePendulumAnglesFromMasses(masses[0], masses[1]);
	double baseTheta1 = baseAngles.front();
	double baseTheta2 = baseAngles.back();

//	std::cout << "Printing Base pendulum: " << std::endl << std::endl;
//	std::cout << "Angles: theta1: " << baseAngles.front() << " theta2: " << baseAngles.back() << std::endl << std::endl;
//
//	std::cout << "Printing renormalized mass properties" << std::endl << std::endl;

//	std::cout << "Mass 1:" << std::endl;
//	std::cout << "Pos: x: " << masses[0].pos.x << " y: " << masses[0].pos.y << std::endl;
//	std::cout << "Vel: x: " << masses[0].linear_vel.x << " y: " << masses[0].linear_vel.y << std::endl;
//	std::cout << "Forces: x: " << masses[0].force_ext.x << " y: " << masses[0].force_ext.y << std::endl << std::endl;
//
//	std::cout << "Mass 2:" << std::endl;
//	std::cout << "Pos: x: " << masses[1].pos.x << " y: " << masses[1].pos.y << std::endl;
//	std::cout << "Vel: x: " << masses[1].linear_vel.x << " y: " << masses[1].linear_vel.y << std::endl;
//	std::cout << "Forces: x: " << masses[1].force_ext.x << " y: " << masses[1].force_ext.y << std::endl << std::endl;

	// Loop over each perturbed pendulum, skip first two masses
//	std::cout << "Normalizing perturbed Pendulums" << std::endl << std::endl;
	int numPerturbedMasses = 2 * numPerturbations;
	for (int i = 0; i < numPerturbations; ++i) {

		// Mass indexes
		int mass1Ind = 2 * i + 2;
		int mass2Ind = mass1Ind + 1;

		// Calculate angles from masses
		std::vector<double> perturbedAngles = DoublePendulumAnglesFromMasses(masses[mass1Ind], masses[mass2Ind]);

//		std::cout << "Perturbed pendulum: " << i << " has masses: " << mass1Ind << " and " << mass2Ind << std::endl;

		// Find the norm of the tangent vector with respect to the basePendulum
		double perturbedTheta1 = perturbedAngles.front();
		double perturbedTheta2 = perturbedAngles.back();

//		std::cout << "Theta1: " << perturbedTheta1 << " Theta2: " << perturbedTheta2 << std::endl;

		Vector tangent = Vector(perturbedTheta1 - baseTheta1, perturbedTheta2 - baseTheta2);
		double tangentNorm = vector_magnitude(tangent);

//		std::cout << "Tangent vector: deltaTheta1: " << tangent.x << " deltaTheta2: " << tangent.y << std::endl;
//		std::cout << "Tangent vector magnitude: " << tangentNorm << std::endl;

		// Calculate sum and store it in lyapunovSums array at correct bucket
		double lyapunovSum = log(tangentNorm / initialSeparation);
		lyapunovRunningSums[i] += lyapunovSum;

		// Renormalize the tangent vector
		tangent.x = tangent.x * initialSeparation / tangentNorm;
		tangent.y = tangent.y * initialSeparation / tangentNorm;

//		std::cout << "Unit vector of tangent vector is: theta1: " << tangent.x / tangentNorm << " theta2: " << tangent.y / tangentNorm << std::endl;
//		std::cout << "Unit vector scaled is: theta1: " << tangent.x / tangentNorm * initialSeparation << " theta2: " << tangent.y / tangentNorm * initialSeparation << std::endl;

		// Apply tangent vector to masses to renormalize angle space position of masses of perturbed pendulum
		double newAngles[2] = {baseTheta1 + tangent.x, baseTheta2 + tangent.y};

//		std::cout << "Renormalized vector: theta1: " << newAngles[0] << " theta2: " << newAngles[1] << std::endl;
		initializeDoublePendulumMassPositions(masses[mass1Ind], masses[mass2Ind], newAngles);

		std::vector<double> renormalizedAngles = DoublePendulumAnglesFromMasses(masses[mass1Ind], masses[mass2Ind]);

		// DEBUGGING
//		std::cout << "Re normalized angles: theta1: " << renormalizedAngles.front() << " theta2: " << renormalizedAngles.back() << std::endl;
//		std::cout << "Printing renormalized mass properties" << std::endl << std::endl;
//		std::cout << "Mass 1:" << std::endl;
//		std::cout << "Pos: x: " << masses[mass1Ind].pos.x << " y: " << masses[mass1Ind].pos.y << std::endl;
//		std::cout << "Vel: x: " << masses[mass1Ind].linear_vel.x << " y: " << masses[mass1Ind].linear_vel.y << std::endl;
//		std::cout << "Forces: x: " << masses[mass1Ind].force_ext.x << " y: " << masses[mass1Ind].force_ext.y << std::endl << std::endl;
//
//		std::cout << "Mass 2:" << std::endl;
//		std::cout << "Pos: x: " << masses[mass2Ind].pos.x << " y: " << masses[mass2Ind].pos.y << std::endl;
//		std::cout << "Vel: x: " << masses[mass2Ind].linear_vel.x << " y: " << masses[mass2Ind].linear_vel.y << std::endl;
//		std::cout << "Forces: x: " << masses[mass2Ind].force_ext.x << " y: " << masses[mass2Ind].force_ext.y << std::endl << std::endl;
//		std::cout << std::endl << std::endl << std::endl;
	}
}

// Evaluate largest lyapunov exponent
double largestLyapunovResult(double *lyapunovRunningSums, int numPerturbations, const int totalSteps, const double dt) {

	// Find maximum running sum
	double* maxLyapunovSum = std::max_element(lyapunovRunningSums, lyapunovRunningSums + numPerturbations);

	// Return maximum lyapunov exponent
	double timeElapsed = (double)totalSteps * dt;
	double maxLyapunovExponent = (*maxLyapunovSum) / timeElapsed;

	return maxLyapunovExponent;
}


// Simulation function
double singleDoublePendulumChaosSim(double baseAngles[2], int numPerturbations, double perturbationRadius, double dt, int totalTime, int framesPerSecond, int lyapunovFrequency, bool data) {

	// Create the base Double Pendulum using default values
	DoublePendulum basePendulum(baseAngles);

	// Define system characteristics
	int systemSize = numPerturbations + 1;
	int numMasses = 2 * systemSize;
	int numBar1s = systemSize;
	int numBar2s = systemSize;

	// Generate the component arrays
	std::tuple<Circular_Rigid_Body *, Rigid_Bar_1 *, Rigid_Bar_2 *> simData = CreateComponentArrays(basePendulum, baseAngles, perturbationRadius, numPerturbations);
	Circular_Rigid_Body *mass_list = std::get<0>(simData);
	Rigid_Bar_1 *bar1s = std::get<1>(simData);
	Rigid_Bar_2 *bar2s = std::get<2>(simData);

	// Set initial conditions for bars
	for (int i = 0; i < systemSize; ++i) {
		bar1s[i].determine_initial_point(mass_list);
		bar2s[i].determine_initial_points(mass_list);
	}

	// Total Simulation Timesteps
	int totalSteps = totalTime / dt;

	// Position data collection: Store the x,y positions of each masses
	// Pivot point is assumed <0,0>
	int stepsPerSec = std::ceil(1 / dt);
	int ITERATIONS = std::ceil(stepsPerSec / framesPerSecond);

	// Collect data every ITERATIONS step
	int dataLen = ((totalSteps + ITERATIONS - 1) / ITERATIONS);

	// Data for masses of each pendulum
	Vector **massData = (Vector **)malloc(dataLen * sizeof(Vector *));
	for (int i = 0; i < dataLen; ++i) {
		massData[i] = (Vector *)malloc(numMasses * sizeof(Vector));
	}

	// Simulation loop
	std::cout << "Total_Steps: " << totalSteps << std::endl;
	int loopCount = 0;
	// Data bucket count
	int dc = 0;

	// Initial energy
	double initialEnergy = total_energy(mass_list, numMasses);
//	std::cout << "Initial energy: " << initialEnergy << std::endl;

	// Lyapunov running sum array
	double *lyapunovSums = new double[systemSize]();

	// Data taking loop
	// Collects data in array and outputs to csv at end
	if (data) {
		while (loopCount < totalSteps)
		{
			// Report energy conservation
			if (loopCount % 300 == 0) {
				std::cout << "Loop: " << loopCount << std::endl;
				double currentEnergy = total_energy(mass_list, numMasses);
				std::cout << "Current energy: " << currentEnergy << std::endl;
				std::cout << "Percent error: " << ((currentEnergy - initialEnergy) / initialEnergy) * 100.0 << std::endl << std::endl;
			}

			// Store pendulum data
			if (loopCount % ITERATIONS == 0)
			{
				// Masses
				for (int i = 0; i < numMasses; ++i)
				{
					massData[dc][i] = {mass_list[i].pos.x, mass_list[i].pos.y};
				}

				// Increment data bucket count
				dc++;
			}

			// Perform lyapunov procedure over masses
			if (loopCount % lyapunovFrequency == 0) {
				// TIMESTEP LYAPUNOV STUFF
				singleDoublePendulumTimestepLyapunov(mass_list, lyapunovSums, numPerturbations, perturbationRadius);
			}

			// rk4 method
			rk4(mass_list, numMasses, bar1s, numBar1s, bar2s, numBar2s, dt);

			// DEBUGGING
//			if (loopCount % 3) {
//				std::cout << "Loop: " << loopCount << std::endl << std::endl;
//
//				// Base pendulum debugging: ALL GOOD
//				std::cout << "Printing Base Pendulum properties" << std::endl;
//				std::vector<double> angs = DoublePendulumAnglesFromMasses(mass_list[0], mass_list[1]);
//				std::cout << "Angles: theta1: " << angs.front() << " theta2: " << angs.back() << std::endl;
////				std::cout << "Rigid Bar constraints" << std::endl;
////				std::cout << "Bar1: " << bar1s[0].constraint(mass_list) << std::endl;
////				std::cout << "Bar2: " << bar2s[0].constraint(mass_list) << std::endl << std::endl;
////
////				std::cout << "Printing base pendulum masses" << std::endl << std::endl;
////
////				std::cout << "Mass: " << 0 << std::endl;
////				std::cout << "Pos: x: " << mass_list[0].pos.x << " y: " << mass_list[0].pos.y << std::endl;
////				std::cout << "Vel: x: " << mass_list[0].linear_vel.x << " y: " << mass_list[0].linear_vel.y << std::endl;
////				std::cout << "Forces: x: " << mass_list[0].force_ext.x << " y: " << mass_list[0].force_ext.y << std::endl << std::endl;
////
////				std::cout << "Mass: " << 1 << std::endl;
////				std::cout << "Pos: x: " << mass_list[1].pos.x << " y: " << mass_list[1].pos.y << std::endl;
////				std::cout << "Vel: x: " << mass_list[1].linear_vel.x << " y: " << mass_list[1].linear_vel.y << std::endl;
////				std::cout << "Forces: x: " << mass_list[1].force_ext.x << " y: " << mass_list[1].force_ext.y << std::endl << std::endl;
//
//				std::cout << "Perturbed pendulum: " << 0 << std::endl;
//				angs = DoublePendulumAnglesFromMasses(mass_list[2], mass_list[3]);
//				std::cout << "Angles: theta1: " << angs.front() << " theta2: " << angs.back() << std::endl << std::endl;
//
////				std::cout << "Rigid Bar constraints" << std::endl;
////				std::cout << "Bar1: " << bar1s[1].constraint(mass_list) << std::endl;
////				std::cout << "Bar2: " << bar2s[1].constraint(mass_list) << std::endl << std::endl;
////				std::cout << "Printing properties of masses: " << 2 << " and: " << 3 << std::endl << std::endl;
////				std::cout << "Mass: " << 2 << std::endl;
////				std::cout << "Pos: x: " << mass_list[2].pos.x << " y: " << mass_list[2].pos.y << std::endl;
////				std::cout << "Vel: x: " << mass_list[2].linear_vel.x << " y: " << mass_list[2].linear_vel.y << std::endl;
////				std::cout << "Mass: " << 3 << std::endl;
////				std::cout << "Pos: x: " << mass_list[3].pos.x << " y: " << mass_list[3].pos.y << std::endl;
////				std::cout << "Vel: x: " << mass_list[3].linear_vel.x << " y: " << mass_list[3].linear_vel.y << std::endl << std::endl << std::endl;
//
////				// Perturbed pendulum debugging
////				for (int i = 2; i < numMasses - 1; i += 2) {
////					std::cout << "Perturbed pendulum: " << (i - 2) / 2 << std::endl << std::endl;
////					std::cout << "Rigid Bar constraints" << std::endl;
////					std::cout << "Bar1: " << bar1s[0].constraint(mass_list) << std::endl;
////					std::cout << "Bar2: " << bar2s[0].constraint(mass_list) << std::endl << std::endl;
////					std::cout << "Printing properties of masses: " << i << " and: " << i + 1 << std::endl << std::endl;
////					std::cout << "Mass: " << i << std::endl;
////					std::cout << "Pos: x: " << mass_list[i].pos.x << " y: " << mass_list[i].pos.y << std::endl;
////					std::cout << "Vel: x: " << mass_list[i].linear_vel.x << " y: " << mass_list[i].linear_vel.y << std::endl;
////
////					std::cout << "Mass: " << i + 1 << std::endl;
////					std::cout << "Pos: x: " << mass_list[i + 1].pos.x << " y: " << mass_list[i + 1].pos.y << std::endl;
////					std::cout << "Vel: x: " << mass_list[i + 1].linear_vel.x << " y: " << mass_list[i + 1].linear_vel.y << std::endl;
////				}
//			}

			// Increment the loop count
			loopCount++;
		}

		// Write angle and position data to file
		DoublePendulumPerturbationsToAnglesCSV(massData, dataLen, numMasses);
		DoublePendulumPerturbationsToCSV(massData, dataLen, numMasses);
	}
	// Non data taking loop
	else {
		while (loopCount < totalSteps)
		{
			// Report energy conservation
			if (loopCount % 300 == 0) {
				std::cout << "Loop: " << loopCount << std::endl;
				double currentEnergy = total_energy(mass_list, numMasses);
				std::cout << "Current energy: " << currentEnergy << std::endl;
				std::cout << "Percent error: " << ((currentEnergy - initialEnergy) / initialEnergy) * 100.0 << std::endl << std::endl;
			}

			// Perform lyapunov procedure over masses
			if (loopCount % lyapunovFrequency == 0) {
				// TIMESTEP LYAPUNOV STUFF
				singleDoublePendulumTimestepLyapunov(mass_list, lyapunovSums, numPerturbations, perturbationRadius);
			}
			// rk4 method
			rk4(mass_list, numMasses, bar1s, numBar1s, bar2s, numBar2s, dt);
			// Increment the loop count
			loopCount++;
		}
	}

	// Evaluate lyapunov exponents
//	double maxLyapunov = largestLyapunovResult(lyapunovSums, numPerturbations, totalSteps, dt);

	// Success
	std::cout << "Done" << std::endl;

	// Free the memory
	free(mass_list);
	free(bar1s);
	free(bar2s);
	delete[] lyapunovSums;

	// Return largest lyapunov
	return 0;
//	return maxLyapunov;
}


