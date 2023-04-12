/*
 * DoublePendulumChaos.cpp
 *
 *  Created on: Apr 11, 2023
 *      Author: norms
 */

#include <math.h>
#include <vector>
#include <tuple>
#include <iomanip>
#include <stdlib.h>

#include "../headers/DoublePendulumChaos.h"
#include "../headers/rigid_bodies.h"
#include "../headers/constraint_bodies.h"
#include "../headers/DoublePendulum.h"

// WORKING
// Function that returns flattened array of angles corresponding to basePendulums evenly spaced in a rectangular grid in the angle space.
double *basePendulumAngleGrid(double topLeft[2], double bottomRight[2], int resolution[2]) {

	// Determine the number of points in each axis
    double theta1Dist = bottomRight[0] - topLeft[0];
    double theta2Dist = topLeft[1] - bottomRight[1];

    // Calculate the step size for each axis
    double xStep = theta1Dist / (double)resolution[0];
    double yStep = theta2Dist / (double)resolution[1];

//    std::cout << "Theta1 distance: " << bottomRight[0] - topLeft[0] << std::endl;
//    std::cout << "Theta2 distance: " << topLeft[1] - bottomRight[1] << std::endl;
//    std::cout << "xStep: " << xStep << " yStep: " << yStep << std::endl;

    // Allocate memory for the flattened array
    double *angles = (double *)malloc(2 * resolution[0] * resolution[1] * sizeof(double));

    // Index that increments along the flattened array
    int linearIndex = 0;
    // Fill in the flattened array with angles
    for (int i = 0; i < resolution[0]; i++)
    {
        for (int j = 0; j < resolution[1]; j++)
        {
            // Calculate the angle for the current position on the grid
        	angles[linearIndex++] = topLeft[0] + i * xStep;
        	angles[linearIndex++] = bottomRight[1] + j * yStep;
        }
    }

    // Return the pointer to the flattened array of angles
    return angles;
}


// Testing
// Function that reads over all of the angles, creates giant arrays full of pendulums and perturbations, rigidBar1 array, and rigidBar2 array. Only perturbed in angle space
std::tuple<Circular_Rigid_Body*, Rigid_Bar_1*, Rigid_Bar_2*> bigArrayMaker(double *baseAngles, int numBasePendulums, int numPerturbations, double radius) {
    std::cout << "Starting bigArrayMaker function." << std::endl;

    // Find total number of pendulums
    int numPerturbedPendulums = numBasePendulums * numPerturbations;
    int numPendulums = numBasePendulums + numPerturbedPendulums;

    // System size determines how many to skip
    int systemSize = numPerturbations + 1;
    int massesPerSystem = 2 * systemSize;

    // Initialize big arrays

    // Masses: First half holds all of the base Pendulum masses, second half holds all of the perturbed pendulum masses
    Circular_Rigid_Body *masses = (Circular_Rigid_Body *)malloc(2 * numPendulums * sizeof(Circular_Rigid_Body));

    // Constraints
    Rigid_Bar_1 *bar1s = (Rigid_Bar_1 *)malloc(numPendulums * sizeof(Rigid_Bar_1));
    Rigid_Bar_2 *bar2s = (Rigid_Bar_2 *)malloc(numPendulums * sizeof(Rigid_Bar_2));

    // Fillup array with offsets of each perturbedAngle calculation
    double perturbedAngles[2 * numPerturbations];

    // Preprocessing
    double c = 2 * M_PI / numPerturbations;

    for (int i = 0; i < numPerturbations; ++i) {
        // Preprocessing
        double angle = c * i;
        int ind1 = 2 * i;
        int ind2 = ind1 + 1;

        // Store most of angle
        perturbedAngles[ind1] = radius * cos(angle);
        perturbedAngles[ind2] = radius * sin(angle);
    }

    // Initialization stuff
    Vector pivot = {0, 0};

    // Iterate over all basePendulums

    for (int i = 0; i < numBasePendulums; ++i) {

        // Initializations:
    	int mass1Ind = i * massesPerSystem;
    	int mass2Ind = mass1Ind + 1;
    	int barInd = 2 * systemSize;
    	initializeDefaultDoublePendulumMasses(masses[mass1Ind], masses[mass2Ind], baseAngles[2 * i], baseAngles[2 * i + 1]);
    	bar1s[barInd] = Rigid_Bar_1(pivot, mass1Ind);
    	bar2s[barInd] = Rigid_Bar_2(mass1Ind, mass1Ind + 1);

    	for (int j = 0; j < numPerturbations; ++j) {

    		// Index processing
    		int mass1Ind = (i * massesPerSystem) + 2 + (2 * j);
    		int mass2Ind = mass1Ind + 1;

    		// Calculate perturbed angles from angle offset array and baseAngles
    		int angle1Ind = 2 * j;
    		int angle2Ind = angle1Ind + 1;
    		int baseAngle1Ind = 2 * i;
    		int baseAngle2Ind = baseAngle1Ind + 1;


    		double angle1 = baseAngles[baseAngle1Ind] + perturbedAngles[angle1Ind];
    		double angle2 = baseAngles[baseAngle2Ind] + perturbedAngles[angle2Ind];

    		int barInd = (i * systemSize) + j + 1;

    		initializeDefaultDoublePendulumMasses(masses[mass1Ind], masses[mass2Ind], angle1, angle2);
    		bar1s[barInd] = Rigid_Bar_1(pivot, mass1Ind);
    		bar2s[barInd] = Rigid_Bar_2(mass1Ind, mass2Ind);

    		// DEBUGGING
    		std::cout << "j: " << j << std::endl;
			std::cout << "mass1Ind: " << mass1Ind << std::endl;
			std::cout << "mass2Ind: " << mass2Ind << std::endl;
			std::cout << "angle1Ind: " << angle1Ind << std::endl;
			std::cout << "angle2Ind: " << angle2Ind << std::endl;
			std::cout << "baseAngle1Ind: " << baseAngle1Ind << std::endl;
			std::cout << "baseAngle2Ind: " << baseAngle2Ind << std::endl;
			std::cout << "angle1: " << angle1 << std::endl;
			std::cout << "angle2: " << angle2 << std::endl << std::endl;
    	}
    }

	// Return the tuple of pointers to the arrays
	return std::make_tuple(masses, bar1s, bar2s);
}
// Working function
// Function that reads over all of the angles, creates giant arrays full of pendulums and perturbations, rigidBar1 array, and rigidBar2 array. Only perturbed in angle space
// Function that reads over all of the angles, creates giant arrays full of pendulums and perturbations, rigidBar1 array, and rigidBar2 array. Only perturbed in angle space
std::tuple<Circular_Rigid_Body*, Rigid_Bar_1*, Rigid_Bar_2*> bigArrayMaker_old(double *baseAngles, int numBasePendulums, int numPerturbations, double radius) {
    std::cout << "Starting bigArrayMaker function." << std::endl;

    // Find total number of pendulums
    int numPerturbedPendulums = numBasePendulums * numPerturbations;
    int numPendulums = numBasePendulums + numPerturbedPendulums;

    // Initialize big arrays
    int numMasses = 2 * numPendulums;
    int numBar1s = numPendulums;
    int numBar2s = numPendulums;

    // Masses: First half holds all of the base Pendulum masses, second half holds all of the perturbed pendulum masses
    Circular_Rigid_Body *masses = (Circular_Rigid_Body *)malloc(numMasses * sizeof(Circular_Rigid_Body));

    // Constraints
    Rigid_Bar_1 *bar1s = (Rigid_Bar_1 *)malloc(numBar1s * sizeof(Rigid_Bar_1));
    Rigid_Bar_2 *bar2s = (Rigid_Bar_2 *)malloc(numBar2s * sizeof(Rigid_Bar_2));

    // Fillup array with offsets of each perturbedAngle calculation
    double perturbedAngles[2 * numPerturbations];

    // Preprocessing
    double c = 2 * M_PI / numPerturbations;

    for (int i = 0; i < numPerturbations; ++i) {
        // Preprocessing
        double angle = c * i;
        int ind1 = 2 * i;
        int ind2 = ind1 + 1;

        // Store most of angle
        perturbedAngles[ind1] = radius * cos(angle);
        perturbedAngles[ind2] = radius * sin(angle);
    }



    // Initialization stuff
    Vector pivot = {0, 0};

    // Iterate over all basePendulums
    for (int i = 0; i < numBasePendulums; ++i) {

        // Index preprocessing
        int basePendulumMass1Index = 2 * i;
        int basePendulumMass2Index = basePendulumMass1Index + 1;

        // Initializations:
    	initializeDefaultDoublePendulumMasses(masses[basePendulumMass1Index], masses[basePendulumMass2Index], baseAngles[basePendulumMass1Index], baseAngles[basePendulumMass2Index]);
    	bar1s[i] = Rigid_Bar_1(pivot, basePendulumMass1Index);
    	bar2s[i] = Rigid_Bar_2(basePendulumMass1Index, basePendulumMass2Index);

    }

    // Iterate over all systems of perturbed pendulums
    int perturbedOffset = numBasePendulums;
    int perturbedMassOffset = 2 * numBasePendulums;

    // Iterate over all groups of perturbed pendulums, 1 group per basePendulum
    for (int i = 0; i < numBasePendulums; ++i) {

    	// Iterate over each perturbed pendulum
    	for (int j = 0; j < numPerturbations; ++j) {

    		int perturbedAngle1Ind = 2 * j;
    		int perturbedAngle2Ind = perturbedAngle1Ind + 1;


			int mass1Ind = perturbedMassOffset + 2 * (i * numPerturbations + j);
			int mass2Ind = mass1Ind + 1;

			// References to the baseAngles
			int angle1Ind = 2 * i;
			int angle2Ind = angle1Ind + 1;

			double angle1 = baseAngles[angle1Ind] + perturbedAngles[perturbedAngle1Ind];
			double angle2 = baseAngles[angle2Ind] + perturbedAngles[perturbedAngle2Ind];

			initializeDefaultDoublePendulumMasses(masses[mass1Ind], masses[mass2Ind], angle1, angle2);

			int perturbedBarInd = perturbedOffset + i * numPerturbations + j;
			bar1s[perturbedBarInd] = Rigid_Bar_1(pivot, mass1Ind);
			bar2s[perturbedBarInd] = Rigid_Bar_2(mass1Ind, mass2Ind);
    	}
    }

	// Return the tuple of pointers to the arrays
	return std::make_tuple(masses, bar1s, bar2s);
}




