/*
 * toFile.h
 *
 *  Created on: Feb 21, 2023
 *      Author: norms
 */

#ifndef HEADERS_TOFILE_H_
#define HEADERS_TOFILE_H_

#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <sstream>
#include <vector>
#include "vectors.h"
#include "rigid_bodies.h"

// Writes the pendulum position data to file. The pivot position, and positions of all masses
void pendulumToFile(Vector pivotData, Vector **massData, int dataLen, int num_masses);

// Write the grid of angles to csv
void writeAnglesToCSV(double *angles, int resolution[2]);

// Write initial angles of all pendulums initialized about the grid of angles in the simulation to csv for plotting.
std::vector<std::vector<double>> DoublePendulumAnglesFromMassListAndWriteToCSV(Circular_Rigid_Body *mass_list, int mass_list_size);

#endif /* HEADERS_TOFILE_H_ */
