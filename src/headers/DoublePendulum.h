/*
 * Pendulum.h
 *
 *  Created on: Feb 24, 2023
 *      Author: norms
 */

#ifndef HEADERS_DOUBLEPENDULUM_H_
#define HEADERS_DOUBLEPENDULUM_H_

#include "rigid_bodies.h"
#include "get_state.h"
#include "utils.h"
#include <vector>
#include <fstream>
#include <sstream>

// Class for a double pendulum

class DoublePendulum {
public:
	Circular_Rigid_Body masses[2];
	Rigid_Bar_1 bar1;
	Rigid_Bar_2 bar2;
	double barLen1;
	double barLen2;

	// Constructors
	DoublePendulum();
	DoublePendulum(double angles[2]);
    DoublePendulum(double angles[2], double barLens[2], double massAmounts[2], Vector vels[2], Vector forces[2], int massInds[2]);
    virtual ~DoublePendulum();

    // Initializes the circular rigid body objects in double pendulum based off params provided in constructor
    void initializeMasses(double angles[2], double massAmounts[2], Vector vels[2], Vector forces[2]);
    // Returns vector of perturbed pendulums
    std::vector<DoublePendulum> *createPerturbedPendulums(int numPerturbations, double radius);
    // Returns the angles based on the cartesian positions of the masses
    std::vector<double> getAngles();
    // ... other methods
};

// ---------------------------------- Functions separated from the class that deal with the parts that construct it      -----------------------------------------------------
std::tuple<Circular_Rigid_Body*, Rigid_Bar_1*, Rigid_Bar_2*> getSimulationData(DoublePendulum &basePendulum, std::vector<DoublePendulum> &perturbedPendulums);
void DoublePendulumSim(double baseAngles[2], int numPerturbations, double perturbationRadius, double dt, int totalTime, int framesPerSecond, bool data);

// Create component arrays for a system of one base pendulum and a number of perturbed pendulums
std::tuple<Circular_Rigid_Body*, Rigid_Bar_1*, Rigid_Bar_2*> CreateComponentArrays(DoublePendulum &basePendulum, double *baseAngles, double radius, int numPerturbations);

// Initializes two CircularMassObjects given parameters. Does not have to be linked to a Double Pendulum Class
void initializeDoublePendulumMasses(Circular_Rigid_Body &mass1, Circular_Rigid_Body &mass2, double *angles, double *barLens, double *massAmounts, Vector *vels, Vector *forces);

// Initialize only the double pendulum masses given a set of two angles
void initializeDoublePendulumMassPositions(Circular_Rigid_Body &mass1, Circular_Rigid_Body &mass2, double *angles);

// Initializes a default set of double pendulum masses given a set of two angles
void initializeDefaultDoublePendulumMasses(Circular_Rigid_Body &mass1, Circular_Rigid_Body &mass2, double angle1, double angle2);

// Completely initiales a set of double pendulum masses given all information
std::vector<double> DoublePendulumAnglesFromMasses(Circular_Rigid_Body &mass1, Circular_Rigid_Body &mass2);

// Functions that write to file

// Writes position data of base pendulum and perturbed pendulums to csv
void DoublePendulumPerturbationsToCSV(Vector** massData, int dataLen, int numMasses);

// Writes angle data (in degrees) of base pendulum and perturbed pendulums to csv
void DoublePendulumPerturbationsToAnglesCSV(Vector** massData, int dataLen, int numMasses);
// Writes initial angles of base pendulum and perturbed pendulums to csv
void DoublePendulumPerturbationInitialAnglesToCSV(DoublePendulum &basePendulum, std::vector<DoublePendulum> &perturbedPendulums);
#endif /* HEADERS_DOUBLEPENDULUM_H_ */
