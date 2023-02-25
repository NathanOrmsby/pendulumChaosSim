/*
 * singlePendulumSim.h
 *
 *  Created on: Feb 24, 2023
 *      Author: norms
 */

#ifndef HEADERS_SINGLEPENDULUMSIM_H_
#define HEADERS_SINGLEPENDULUMSIM_H_

#include "rigid_bodies.h"
#include "constraint_bodies.h"

// Simulates a single pendulum of specified customization for a specified amount of time
void singlePendulumSim(int numMasses, Circular_Rigid_Body *masses, int numBar1s, Rigid_Bar_1 *bar1s,  int numBar2s, Rigid_Bar_2 *bar2s, double dt, int totalTime);



#endif /* HEADERS_SINGLEPENDULUMSIM_H_ */
