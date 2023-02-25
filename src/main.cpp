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
#include "headers/get_state.h"
#include "headers/utils.h"
#include "headers/renderer.h"
#include "headers/matrix_stuff.h"
#include "headers/rk4.h"
#include "headers/rigid_bodies.h"
#include "headers/springs.h"
#include "headers/toFile.h"
#include "headers/total_energy.h"

// Main function headers'
#include "headers/singlePendulumSim.h"

// Globals
bool running = true;

int main(void) {

	// ********* INITIAL CONDITIONS *************

	// Where the initial conditions are set.

	// MASSES

	// For double pendulum
	int num_masses = 2;

	Circular_Rigid_Body masses[num_masses];

	masses[0].pos = {7.32, 4.4};
	masses[1].pos = {7.32, 5.4};

	// Fill in the rest of the parameters for the masses
	for (int i = 0; i < num_masses; i++)
	{
		masses[i].linear_vel = {0.0, 0.0};
		masses[i].force_ext = {0.0, 0.0};
		masses[i].mass = 1.0;
	}
	// Constraints lists
	int num_bar1s = 1;
	int num_bar2s = 1;

	Rigid_Bar_1 bar1s[num_bar1s];
	bar1s[0].pivot = {6.32, 3.4};
	bar1s[0].attached_mass = 0;
	bar1s[0].determine_initial_point(masses);


	Rigid_Bar_2 bar2s[num_bar2s];


	bar2s[0].attached_masses[0] = 0;
	bar2s[0].attached_masses[1] = 1;
	bar2s[0].determine_initial_points(masses);

	// Springs lists
	//int num_spring1s = 0;
	int num_spring2s = 0;

	Spring_2 spring2s[num_spring2s];
//	spring2s[0].attached_masses[0] = 1;
//	spring2s[0].attached_masses[1] = 4;

	for (int i = 0; i < num_spring2s; i++)
	{
		spring2s[i].eq_length = 0.75;
		spring2s[i].spring_const = 5.0;
		spring2s[i].color = 0xFFFF00;
	}


	// Try calling our new function
	double dt = 1.0 / 1000.0;
	int totalTime = 10;
	singlePendulumSim(num_masses, masses, num_bar1s, bar1s,  num_bar2s, bar2s, dt, totalTime);

//	// Create the State_Getter class for current state
//	State_Getter current_state;
//	// Number of masses and constraints
//	current_state.num_bodies = num_masses;
//	current_state.num_constraints = num_bar1s + num_bar2s;
//	current_state.num_spring2s = num_spring2s;
//
//	// Timestep
//	double dt = 1.0 / 1000.0;
//
//	// Total run time in seconds
//	int totalTime = 10;
//
//	// Total Steps
//	int totalSteps = totalTime / dt;
//
//	// Data collection
//	// Data array for storage of position data: Store the pivot position, and mass positions
//
//	// Collect data every ITERATIONS step
//	int ITERATIONS = 17;
//	int dataLen = ((totalSteps + ITERATIONS - 1) / ITERATIONS) + 1;
//
//	// Data for pivot point
//	Vector pivotData = {bar1s[0].pivot.x, bar1s[0].pivot.y};
//
//
//	// Data for Masses
//	Vector **massData = (Vector **)malloc(dataLen * sizeof(Vector *));
//
//	// Allocate space for the masses
//	for (int i = 0; i < dataLen; ++i)
//	{
//		massData[i] = (Vector *)malloc(num_masses * sizeof(Vector));
//	}
//
//
//	// Loop
//
//	// Loop count
//	int loopCount = 0;
//	// Data bucket count
//	int dc = 0;
//	while (loopCount < totalSteps)
//	{
//		// Store data
//		if (loopCount % ITERATIONS == 0)
//		{
//			// Masses
//			for (int i = 0; i < num_masses; ++i)
//			{
//				massData[dc][i] = {masses[i].pos.x, masses[i].pos.y};
//			}
//
//			// Increment data bucket count
//			dc++;
//
//			std::cout << "Constraints" << std::endl;
//			std::cout << "c1: " << bar1s[0].constraint(masses) << std::endl;
//			std::cout << "c2: " << bar2s[0].constraint(masses) << std::endl;
//			std::cout << "Total energy: " << total_energy(masses, num_masses, loopCount) << std::endl << std::endl;
//
//		}
//
//		// Numerical integrator
//		// rk4 method
//		rk4(&current_state, masses, bar1s, num_bar1s, bar2s, num_bar2s, spring2s, dt);
//
//		// Increment the loop count
//		loopCount++;
//	}
//
//	// Record last bucket of data for masses
//	for (int i = 0; i < num_masses; ++i)
//	{
//		massData[dataLen - 1][i] = {masses[i].pos.x, masses[i].pos.y};
//	}
//
//
//	// Write data to file
//	pendulumToFile(pivotData, massData, dataLen, num_masses);
//
//	// Success
//	std::cout << "Done" << std::endl;
//	return 0;
}
