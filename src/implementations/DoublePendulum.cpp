/*
 * Pendulum.cpp
 *
 *  Created on: Feb 24, 2023
 *      Author: norms
 */

#include "../headers/Vectors.h"
#include "../headers/rk4.h"
#include "../headers/total_energy.h"

#include <math.h>
#include <vector>
#include <tuple>
#include <iomanip>
#include "../headers/DoublePendulum.h"
#include <iomanip>
#include "../headers/toFile.h"


using namespace std;

// Methods that act on the Pendulum Classes

// Default constructor
DoublePendulum::DoublePendulum() : masses{Circular_Rigid_Body(), Circular_Rigid_Body()}, bar1(), bar2(), barLen1(1.0), barLen2(1.0) {}

DoublePendulum::DoublePendulum(double angles[2]) :
    masses(),
    bar1(),
    bar2(),
    barLen1(1.0),
    barLen2(1.0)
{
    // Default values for mass amounts, velocities, and forces
    double massAmounts[2] = {1, 1};
    Vector vels[2] = {Vector(), Vector()};
    Vector forces[2] = {Vector(), Vector()};

    // Initialize masses and bars
    initializeMasses(angles, massAmounts, vels, forces);

    bar1.attached_mass = 0;
    bar1.pivot.x = 0;
    bar1.pivot.y = 0;
    bar1.determine_initial_point(masses);

    bar2.attached_masses[0] = 0;
    bar2.attached_masses[1] = 1;
    bar2.determine_initial_points(masses);
}

DoublePendulum::DoublePendulum(double angles[2], double barLens[2], double massAmounts[2], Vector vels[2], Vector forces[2], int massInds[2]) {
    barLen1 = barLens[0];
    barLen2 = barLens[1];

    // Initialize masses
    initializeMasses(angles, massAmounts, vels, forces);

    // Initialize bars
    Vector pivot = {0, 0};

    bar1 = Rigid_Bar_1(pivot, massInds[0]);
    bar1.determine_initial_point(masses);

    bar2 = Rigid_Bar_2(massInds[0], massInds[1]);
    bar2.determine_initial_points(masses);
}

DoublePendulum::~DoublePendulum() {
}

void DoublePendulum::initializeMasses(double angles[2], double massAmounts[2], Vector vels[2], Vector forces[2]) {

    // Find cartesian coordinates for all masses

	// DEBUGGING
    // std::cout << "Angles: " << angles[0] << " and " << angles[1] << std::endl;
	// std::cout << radiansToDegrees(angles[0] - M_PI) << std::endl;

    // Set masses for Circular Rigid Body classes
    masses[0].mass = massAmounts[0];
    masses[1].mass = massAmounts[1];



    // First mass
    // Calculate position
    masses[0].pos = Vector(barLen1 * sin(angles[0]), barLen1 * cos(angles[0] - M_PI));

    // Offset origin
    Vector offset = {masses[0].pos.x, masses[0].pos.y};

    // Second mass

    // Calculate position
    masses[1].pos = Vector(barLen2 * sin(angles[1]) + offset.x, cos(angles[1] - M_PI) + offset.y);

    // Velocities
    masses[0].linear_vel = Vector(vels[0].x, vels[0].y);
    masses[1].linear_vel = Vector(vels[1].x, vels[1].y);


    // Forces
    masses[0].force_ext = Vector(forces[0].x, forces[0].y);
    masses[1].force_ext = Vector(forces[1].x, forces[1].y);
}

std::vector<double> DoublePendulum::getAngles() {

	double x1 = masses[0].pos.x;
	double y1 = masses[0].pos.y;
	double x2 = masses[1].pos.x;
	double y2 = masses[1].pos.y;
	// Axis vector, negative y
	Vector axis = {0, -1};

	// Angle 1 from axis to v1
	Vector v1 = {x1, y1};
	double angle1 = atan2(v1.y, v1.x) - atan2(axis.y, axis.x);

	// Angle 2: From axis to v2
	Vector v2 = {x2 - x1, y2 - y1};
	double angle2 = atan2(v2.y, v2.x) - atan2(axis.y, axis.x);

	// Normalize angles between -pi and pi
	if (angle1 > M_PI) angle1 -= 2 * M_PI;
	else if (angle1 <= -M_PI) angle1 += 2 * M_PI;
	if (angle2 > M_PI) angle2 -= 2 * M_PI;
	else if (angle2 <= -M_PI) angle2 += 2 * M_PI;

	return std::vector<double>{angle1, angle2};

}

std::vector<DoublePendulum> *DoublePendulum::createPerturbedPendulums(int numPerturbations, double radius) {

	// Create a vector of pendulum objects perturbed about the base pendulum in a circle in the angle space theta1,theta2
    std::vector<DoublePendulum> *perturbedPendulums = new std::vector<DoublePendulum>;

    std::vector<double> angles = getAngles();



    double *baseAngles = new double[2]{angles.front(), angles.back()};
    std::cout << "Base angles are: theta1: " << baseAngles[0] << " theta2: " << baseAngles[1] << std::endl;
    double *baseMassAmounts = new double[2]{masses[0].mass, masses[1].mass};
    Vector *baseVels = new Vector[2]{Vector(masses[0].linear_vel.x, masses[0].linear_vel.y), Vector(masses[1].linear_vel.x, masses[1].linear_vel.y)};

//    Vector *baseVels = new Vector[2];
//    baseVels[0] = Vector(masses[0].linear_vel.x, masses[0].linear_vel.y);
//    baseVels[1] = Vector(masses[1].linear_vel.x, masses[0].linear_vel.y);
    Vector *baseForces = new Vector[2]{Vector(masses[0].force_ext.x, masses[0].force_ext.y), Vector(masses[1].force_ext.x, masses[1].force_ext.y)};
    double *baseBarLens = new double[2]{barLen1, barLen2};

    // Dynamically allocated variables.
    double *perturbedAngles = new double[2];
    int *massInds = new int[2];

    // Base pendulum is inserted at start, and contains masses 0 and 1, add an offset of 2 for perturbed pendulum masses
    for (int i = 0; i < numPerturbations; ++i) {
        double angle = 2 * M_PI * i / numPerturbations;
        perturbedAngles[0] = baseAngles[0] + radius * cos(angle);
        perturbedAngles[1] = baseAngles[1] + radius * sin(angle);
        massInds[0] = 2 * i + 2;
        massInds[1] = 2 * i + 3;
        std::cout << "Perturbation: " << i << " Theta1: " << perturbedAngles[0] << " Theta2: " << perturbedAngles[1] << std::endl;
        DoublePendulum perturbedPendulum(perturbedAngles, baseBarLens, baseMassAmounts, baseVels, baseForces, massInds);
        perturbedPendulums->push_back(perturbedPendulum);
    }

    delete[] baseVels;
    delete[] baseBarLens;
    delete[] perturbedAngles;
    delete[] massInds;

    return perturbedPendulums;
}

std::tuple<Circular_Rigid_Body*, Rigid_Bar_1*, Rigid_Bar_2*> getSimulationData(DoublePendulum &basePendulum, std::vector<DoublePendulum> &perturbedPendulums) {

	// Insert the base pendulum at the front of the vector
    perturbedPendulums.insert(perturbedPendulums.begin(), basePendulum);

    // Calculate the necessary length of each array
    int numPendulums = perturbedPendulums.size();
    int numMasses = 2 * numPendulums;
    int numBar1s = numPendulums;
    int numBar2s = numPendulums;

    // Allocate memory for arrays
    Circular_Rigid_Body *mass_list = new Circular_Rigid_Body[numMasses];
    Rigid_Bar_1 *bar1s = new Rigid_Bar_1[numBar1s];
    Rigid_Bar_2 *bar2s = new Rigid_Bar_2[numBar2s];

    // Iterate over the entire vector
    for (int i = 0; i < numPendulums; ++i) {
        DoublePendulum &pendulum = perturbedPendulums[i];

        // Fill up the mass_list array
        mass_list[2 * i] = pendulum.masses[0];
        mass_list[2 * i + 1] = pendulum.masses[1];

        // Fill up the bar1s array
        bar1s[i] = pendulum.bar1;

        // Fill up the bar2s array
        bar2s[i] = pendulum.bar2;
    }

    // Return the tuple of pointers to arrays
    return std::make_tuple(mass_list, bar1s, bar2s);
}

void DoublePendulumSim(double baseAngles[2], int numPerturbations, double perturbationRadius, double dt, int totalTime, int framesPerSecond, bool data)
{
    // Create the base Double Pendulum using default values
    DoublePendulum basePendulum(baseAngles);

//    // DEBUGGING
//    cout << "Printing base pendulum stuff" << endl;
//	cout << "Mass 0: x: " << basePendulum.masses[0].pos.x << " y: " << basePendulum.masses[0].pos.y << endl;
//	cout << "Mass 1: x: " << basePendulum.masses[1].pos.x << " y: " << basePendulum.masses[1].pos.y << endl;
//	std::vector<double> calculatedAngles = basePendulum.getAngles();
//	cout << "Printing angles: theta1: " << radiansToDegrees(calculatedAngles.front()) << " theta2: " << radiansToDegrees(calculatedAngles.back()) << endl;

    // Method that creates redundant pendulum classes
//    // Create the perturbed pendulums
//    std::vector<DoublePendulum> *perturbedPendulums = basePendulum.createPerturbedPendulums(numPerturbations, perturbationRadius);
//    // DEBUGGING: Write angle1,angle2 for each pendulum to csv
////    DoublePendulumPerturbationInitialAnglesToCSV(basePendulum, *perturbedPendulums);
//
//
//    // Generate the tuple of pointers
//    int numMasses = 2 * (numPerturbations + 1);
//    int numBar1s = numPerturbations + 1;
//    int numBar2s = numPerturbations + 1;
//    std::tuple<Circular_Rigid_Body *, Rigid_Bar_1 *, Rigid_Bar_2 *> simulationData = getSimulationData(basePendulum, *perturbedPendulums);
//    Circular_Rigid_Body *mass_list = std::get<0>(simulationData);
//    Rigid_Bar_1 *bar1s = std::get<1>(simulationData);
//    Rigid_Bar_2 *bar2s = std::get<2>(simulationData);


    // Method that doesn't create redundant pendulums
	// Generate the tuple of pointers
	int numMasses = 2 * (numPerturbations + 1);
	int numBar1s = numPerturbations + 1;
	int numBar2s = numPerturbations + 1;
	std::tuple<Circular_Rigid_Body *, Rigid_Bar_1 *, Rigid_Bar_2 *> simData = CreateComponentArrays(basePendulum, baseAngles, perturbationRadius, numPerturbations);
	Circular_Rigid_Body *mass_list = std::get<0>(simData);
	Rigid_Bar_1 *bar1s = std::get<1>(simData);
	Rigid_Bar_2 *bar2s = std::get<2>(simData);

	std::vector<std::vector<double>> angs = DoublePendulumAnglesFromMassListAndWriteToCSV(mass_list, numMasses);

    // Print mass attributes
//    std::cout << "Initial conditions for each pendulum:" << std::endl << std::endl;
//    for (int i = 0; i < numMasses; ++i)
//    {
//        int pendulum_num = i / 2;
//        std::cout << std::endl << "Pendulum " << pendulum_num << ":" << std::endl;
//        std::cout << "Mass " << i << ": " << std::endl
//            << "Pos: (" << mass_list[i].pos.x << ", " << mass_list[i].pos.y << ")" << std::endl
//            << "Vel: (" << mass_list[i].linear_vel.x << ", " << mass_list[i].linear_vel.y << ")" << std::endl
//            << "Force: (" << mass_list[i].force_ext.x << ", " << mass_list[i].force_ext.y << ")" << std::endl
//            << "Mass: " << mass_list[i].mass << std::endl << std::endl;
//    }
//
//    // Print bar1 attributes
//    std::cout << std::endl << "Rigid_Bar_1 attributes for each pendulum:" << std::endl << std::endl;
//    for (int i = 0; i < numBar1s; ++i)
//    {
//        std::cout << "Pendulum " << i << ":" << std::endl
//            << "Rigid_Bar_1 " << i << ": " << std::endl
//            << "Pivot: (" << bar1s[i].pivot.x << ", " << bar1s[i].pivot.y << ")" << std::endl
//            << "Attached Mass: " << bar1s[i].attached_mass << std::endl << std::endl;
//    }
//
//    // Print bar2 attributes
//    std::cout << std::endl << "Rigid_Bar_2 attributes for each pendulum:" << std::endl << std::endl;
//    for (int i = 0; i < numBar2s; ++i)
//    {
//        std::cout << "Pendulum " << i << ":" << std::endl
//            << "Rigid_Bar_2 " << i << ": " << std::endl
//            << "Initial Points: (" << bar2s[i].initial_point0.x << ", " << bar2s[i].initial_point0.y << ") "
//            << "and (" << bar2s[i].initial_point1.x << ", " << bar2s[i].initial_point1.y << ")" << std::endl
//            << "Attached Masses: " << bar2s[i].attached_masses[0] << " and " << bar2s[i].attached_masses[1] << std::endl << std::endl;
//    }

    // Total Simulation Timesteps
    int totalSteps = totalTime / dt;

    // Data collection
	// Data array for storage of position data: Store the pivot position, and mass positions
    // Take data every iterationsSkipped timesteps for 60 fps
    int stepsPerSec = std::ceil(1 / dt);

    int ITERATIONS = std::ceil(stepsPerSec / framesPerSecond);

    std::cout << "stepsPerSec: " << stepsPerSec << std::endl;
    std::cout << "ITERATIONS: " << ITERATIONS << std::endl;

	// Collect data every ITERATIONS step
	int dataLen = ((totalSteps + ITERATIONS - 1) / ITERATIONS) + 1;


	// Data for masses of each pendulum
	Vector **massData = (Vector **)malloc(dataLen * sizeof(Vector *));

	// Pivot point is <0,0>

	if (data) {
		// Allocate space for the masses
		for (int i = 0; i < dataLen; ++i) {
			massData[i] = (Vector *)malloc(numMasses * sizeof(Vector));
		}
	}

	// Loop
	std::cout << "Total_Steps: " << totalSteps << std::endl;
	// Loop count
	int loopCount = 0;
	// Data bucket count
	int dc = 0;

	// Initial energy
	double initialEnergy = total_energy(mass_list, numMasses);
	std::cout << "Initial energy: " << initialEnergy << std::endl;

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


			// rk4 method
			rk4(mass_list, numMasses, bar1s, numBar1s, bar2s, numBar2s, dt);
			// Increment the loop count
			loopCount++;
		}



		// Write data to file
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

			// rk4 method
			rk4(mass_list, numMasses, bar1s, numBar1s, bar2s, numBar2s, dt);
			// Increment the loop count
			loopCount++;
		}
	}

    // Success
    std::cout << "Done" << std::endl;

    // Free the memory
    free(mass_list);
    free(bar1s);
    free(bar2s);
    // delete[] perturbedPendulums;
}

void DoublePendulumPerturbationsToCSV(Vector** massData, int dataLen, int numMasses) {

    // Open the output file
	// fpath:
	std::string fpath = "D:\\pendulumSimulator\\pendulumData\\";
	std::string fileName = "doublePendulumPerturbations.csv";
    std::ofstream outputFile(fpath + fileName);

    // Write the header row
    for (int i = 0; i < numMasses; ++i) {
        outputFile << "m" << i << "x,m" << i << "y";
        if (i < numMasses - 1) {
            outputFile << ",";
        }
    }
    outputFile << std::endl;

    // Write the data rows
    for (int i = 0; i < dataLen; ++i) {
        for (int j = 0; j < numMasses; ++j) {
            outputFile << massData[i][j].x << "," << massData[i][j].y;
            if (j < numMasses - 1) {
                outputFile << ",";
            }
        }
        outputFile << std::endl;
    }

    // Close the output file
    outputFile.close();
}


void DoublePendulumPerturbationsToAnglesCSV(Vector** massData, int dataLen, int numMasses) {
    // Open the output file
    std::string fpath = "D:\\pendulumSimulator\\pendulumData\\";
    std::string fileName = "doublePendulumPerturbationAngles.csv";
    std::ofstream outputFile(fpath + fileName);

    // Write the header row
    for (int i = 0; i < numMasses; ++i) {
        outputFile << "angle" << i;
        if (i < numMasses - 1) {
            outputFile << ",";
        }
    }
    outputFile << std::endl;

    // Write the data rows
    for (int i = 0; i < dataLen; ++i) {

    	// Write all the data in one row
        for (int j = 0; j < numMasses; j += 2) {
            // Convert mass positions to angles
            Circular_Rigid_Body mass1 = Circular_Rigid_Body(massData[i][j]);
            Circular_Rigid_Body mass2 = Circular_Rigid_Body(massData[i][j + 1]);
            std::vector<double> massAngles = DoublePendulumAnglesFromMasses(mass1, mass2);

            // Write angles to file with precision set to 15 decimal places
            outputFile << std::setprecision(15) << radiansToDegrees(massAngles.front()) << "," << radiansToDegrees(massAngles.back());
            if (j < numMasses - 2) {
                outputFile << ",";
            }
        }
        outputFile << std::endl;
    }

    // Close the output file
    outputFile.close();
}





void DoublePendulumPerturbationInitialAnglesToCSV(DoublePendulum &basePendulum, std::vector<DoublePendulum> &perturbedPendulums) {
    std::vector<DoublePendulum> allPendulums = {basePendulum};
    allPendulums.insert(allPendulums.end(), perturbedPendulums.begin(), perturbedPendulums.end());

    std::string fpath = "D:\\pendulumSimulator\\pendulumData\\";
    std::string fileName = "initialPerturbedAngles.csv";
    std::ofstream outputFile(fpath + fileName);
    if (!outputFile.is_open()) {
        std::cerr << "Unable to open file " << fileName << std::endl;
        return;
    }

    // Write the header row
    for (size_t i = 0; i < allPendulums.size(); ++i) {
        outputFile << "pendulum" << i << "Theta1,pendulum" << i << "Theta2";
        if (i < allPendulums.size() - 1) {
            outputFile << ",";
        }
    }
    outputFile << std::endl;

    // Write the data row
    for (size_t i = 0; i < allPendulums.size(); ++i) {
        std::vector<double> angles = allPendulums[i].getAngles();
        outputFile << std::fixed << std::setprecision(6) << radiansToDegrees(angles[0]) << "," << radiansToDegrees(angles[1]);
        if (i < allPendulums.size() - 1) {
            outputFile << ",";
        }
    }
    outputFile << std::endl;

    outputFile.close();
}

// --------------------------------------------- Functions not requiring the class ---------------------------------------------

// Initializes two CircularMassObjects given parameters. Does not have to be linked to a Double Pendulum Class
void initializeDoublePendulumMasses(Circular_Rigid_Body &mass1, Circular_Rigid_Body &mass2, double *angles, double *barLens, double *massAmounts, Vector *vels, Vector *forces) {

    // Find cartesian coordinates for all masses

	// DEBUGGING
    // std::cout << "Angles: " << angles[0] << " and " << angles[1] << std::endl;
	// std::cout << radiansToDegrees(angles[0] - M_PI) << std::endl;

    // Set masses for Circular Rigid Body classes
    mass1.mass = massAmounts[0];
    mass2.mass = massAmounts[1];

    // First mass
    // Calculate position
    mass1.pos = Vector(barLens[0] * sin(angles[0]), barLens[0] * cos(angles[0] - M_PI));

    // Second mass
    // Calculate position
    mass2.pos = Vector(barLens[0] * sin(angles[1]) + mass1.pos.x, cos(angles[1] - M_PI) + mass1.pos.y);

    // Velocities
    mass1.linear_vel = Vector(vels[0].x, vels[0].y);
    mass2.linear_vel = Vector(vels[1].x, vels[1].y);

    // Forces
    mass1.force_ext = Vector(forces[0].x, forces[0].y);
    mass2.force_ext = Vector(forces[1].x, forces[1].y);
}

// This function only initializes the cartesian positions of the masses based on the angles provided. Assumes bar lengths of 1
void initializeDoublePendulumMassPositions(Circular_Rigid_Body &mass1, Circular_Rigid_Body &mass2, double *angles) {

    // Find cartesian coordinates for all masses
    // First mass
    // Calculate position
    mass1.pos = Vector(sin(angles[0]), cos(angles[0] - M_PI));

    // Second mass
    // Calculate position
    mass2.pos = Vector(sin(angles[1]) + mass1.pos.x, cos(angles[1] - M_PI) + mass1.pos.y);
}

// Initializes two CircularMassObjects given parameters. Does not have to be linked to a Double Pendulum Class
void initializeDefaultDoublePendulumMasses(Circular_Rigid_Body &mass1, Circular_Rigid_Body &mass2, double angle1, double angle2) {

	mass1.mass = 1.0;
    mass2.mass = 1.0;

    // Mass positions
    mass1.pos = Vector(sin(angle1), cos(angle1 - M_PI));
    mass2.pos = Vector(sin(angle2) + mass1.pos.x, cos(angle2 - M_PI) + mass1.pos.y);

    // Velocities
    mass1.linear_vel = Vector(0, 0);
    mass2.linear_vel = Vector(0, 0);

    // Forces
    mass1.force_ext = Vector(0, 0);
    mass2.force_ext = Vector(0, 0);
}

std::tuple<Circular_Rigid_Body*, Rigid_Bar_1*, Rigid_Bar_2*> CreateComponentArrays(DoublePendulum &basePendulum, double *baseAngles, double radius, int numPerturbations) {

    // Calculate the necessary length of each array
    int numPendulums = numPerturbations + 1;
    int numMasses = 2 * numPendulums;
    int numBar1s = numPendulums;
    int numBar2s = numPendulums;

    // Allocate memory for arrays
    Circular_Rigid_Body *mass_list = (Circular_Rigid_Body *)malloc(2 * numPendulums * sizeof(Circular_Rigid_Body));
    Rigid_Bar_1 *bar1s = (Rigid_Bar_1 *)malloc(numPendulums * sizeof(Rigid_Bar_1));
    Rigid_Bar_2 *bar2s = (Rigid_Bar_2 *)malloc(numPendulums * sizeof(Rigid_Bar_2));

    // Find angles for each perturbed pendulum


    // Allocate space
    double **perturbedAngles = (double **)malloc(numPerturbations * sizeof(double *));
    for (int i = 0; i < numPerturbations; ++i) {
    	perturbedAngles[i] = (double *)malloc(2 * sizeof(double));
    }

    for (int i = 0; i < numPerturbations; ++i) {
            double angle = 2 * M_PI * i / numPerturbations;
            perturbedAngles[i][0] = baseAngles[0] + radius * cos(angle);
            perturbedAngles[i][1] = baseAngles[1] + radius * sin(angle);
    }

    // Initialize each of the mass objects, and rigid bar objects
    double *barLens = new double[2]{basePendulum.barLen1, basePendulum.barLen2};
    double *massAmounts = new double[2]{1.0, 1.0};
    Vector *vels = new Vector[2]{Vector(), Vector()};
    Vector *forces = new Vector[2]{Vector(), Vector()};
    Vector pivot = {0, 0};

    // Base pendulum
    initializeDoublePendulumMasses(mass_list[0], mass_list[1], baseAngles, barLens, massAmounts, vels, forces);
    bar1s[0] = Rigid_Bar_1(pivot, 0);
    bar2s[0] = Rigid_Bar_2(0, 1);

    // Perturbations
    for (int i = 0; i < numPerturbations; ++i) {
    	int mass1Ind = 2 * i + 2;
    	int mass2Ind = 2 * i + 3;
    	double angles[2] = {perturbedAngles[i][0], perturbedAngles[i][1]};
    	initializeDoublePendulumMasses(mass_list[mass1Ind], mass_list[mass2Ind], angles, barLens, massAmounts, vels, forces);
    	bar1s[i + 1] = Rigid_Bar_1(pivot, mass1Ind);
    	bar2s[i + 1] = Rigid_Bar_2(mass1Ind, mass2Ind);
    }

    // Free stuff
    for (int i = 0; i < numPerturbations; ++i) {
    	free(perturbedAngles[i]);
    }
    free(perturbedAngles);
    delete[] barLens;
    delete[] massAmounts;
    delete[] vels;
    delete[] forces;


    // Return the tuple of pointers to arrays
    return std::make_tuple(mass_list, bar1s, bar2s);
}

std::vector<double> DoublePendulumAnglesFromMasses(Circular_Rigid_Body &mass1, Circular_Rigid_Body &mass2) {
	double x1 = mass1.pos.x;
	double y1 = mass1.pos.y;
	double x2 = mass2.pos.x;
	double y2 = mass2.pos.y;
	// Axis vector, negative y
	Vector axis = {0, -1};

	// Angle 1 from axis to v1
	Vector v1 = {x1, y1};
	double angle1 = atan2(v1.y, v1.x) - atan2(axis.y, axis.x);

	// Angle 2: From axis to v2
	Vector v2 = {x2 - x1, y2 - y1};
	double angle2 = atan2(v2.y, v2.x) - atan2(axis.y, axis.x);

	// Normalize angles between -pi and pi
	if (angle1 > M_PI) angle1 -= 2 * M_PI;
	else if (angle1 <= -M_PI) angle1 += 2 * M_PI;
	if (angle2 > M_PI) angle2 -= 2 * M_PI;
	else if (angle2 <= -M_PI) angle2 += 2 * M_PI;

	return std::vector<double>{angle1, angle2};
}



