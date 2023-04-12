#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <chrono>
#include <time.h>
#include "math.h"
#include <cmath>
#include <stdint.h>



// Custom header files
#include "headers/constraint_bodies.h"
#include "headers/DoublePendulum.h"
#include "headers/get_state.h"
#include "headers/utils.h"
#include "headers/renderer.h"
#include "headers/matrix_stuff.h"
#include "headers/rk4.h"
#include "headers/rigid_bodies.h"
#include "headers/springs.h"
#include "headers/toFile.h"
#include "headers/total_energy.h"
#include "headers/DoublePendulumChaos.h"
#include "headers/singleDoublePendulumChaos.h"

// Main function headers'
// #include "headers/singlePendulumSim.h"

using namespace std;

// Globals
bool running = true;

int main(void) {

//	 DOUBLE PENDULUM SIMULATION STUFF
	double dt = 1.0 / 1800.0;
	int totalTime = 3;
	int framesPerSecond = 900;

	// Specify angles in degrees, converted to radians later
	double theta1 = degreesToRadians(-135);
	double theta2 = degreesToRadians(-135);
	double angles[2] = {theta1, theta2};

	// Perturbations:
	int numPerturbations = 20;
	double radius = 0.001;

	// Specify data taking:
	bool data = true;

	// Testing chaos simulation with single double pendulum and perturbations
	int lyapunovFrequency = 6;
	singleDoublePendulumChaosSim(angles, numPerturbations, radius, dt, totalTime, framesPerSecond, lyapunovFrequency, data);

	// Normal simulation with perturbations
//	DoublePendulumSim(angles, numPerturbations, radius, dt, totalTime, framesPerSecond, data);

//
//	// singlePendulumSim(num_masses, masses, num_bar1s, bar1s,  num_bar2s, bar2s, dt, totalTime);


//	// Unit TESTING FOR CHAOS PLOT FUNCTIONS
//
//	// Grid specifications
//	double topLeft_x = degreesToRadians(-180);
//	double topLeft_y = degreesToRadians(180);
//	double bottomRight_x = degreesToRadians(180);
//	double bottomRight_y = degreesToRadians(-180);
//	int resolution[2] = {20, 20};
//	double topLeft[2] = {topLeft_x, topLeft_y};
//	double bottomRight[2] = {bottomRight_x, bottomRight_y};
//
//	// Perturbations:
//	int numPerturbations = 20;
//	double radius = degreesToRadians(3.0);
//
//	int numBasePendulums = resolution[0] * resolution[1];
//	int numPerturbedPendulums = numBasePendulums * numPerturbations;
//	int numPendulums = numBasePendulums + numPerturbedPendulums;
//	double *baseAngles = basePendulumAngleGrid(topLeft, bottomRight, resolution);
//
//	writeAnglesToCSV(baseAngles, resolution);
//
//	// Testing giant array generation
//	std::tuple<Circular_Rigid_Body*, Rigid_Bar_1*, Rigid_Bar_2*> bigArrays = bigArrayMaker(baseAngles, numBasePendulums, numPerturbations, radius);
//	Circular_Rigid_Body *mass_list = std::get<0>(bigArrays);
//	Rigid_Bar_1 *bar1s = std::get<1>(bigArrays);
//	Rigid_Bar_2 *bar2s = std::get<2>(bigArrays);
//
//	std::cout << "Made it!" << std::endl;
//	// Output to csv
//	std::vector<std::vector<double>> everyAngle = DoublePendulumAnglesFromMassListAndWriteToCSV(mass_list, 2 * numPendulums);
	return 0;
}
