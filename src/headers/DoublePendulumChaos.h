/*
 * DoublePendulumChaos.h
 *
 *  Created on: Apr 11, 2023
 *      Author: norms
 */



#ifndef HEADERS_DOUBLEPENDULUMCHAOS_H_
#define HEADERS_DOUBLEPENDULUMCHAOS_H_

#include "rigid_bodies.h"
#include "constraint_bodies.h"
// Function that returns flatted array of angles evenly spaced in a rectangular grid in the angle space.
double *basePendulumAngleGrid(double topLeft[2], double bottomRight[2], int resolution[2]);
// Function that reads over all of the angles, creates giant arrays full of pendulums and perturbations, rigidBar1 array, and rigidBar2 array. Only perturbed in angle space
std::tuple<Circular_Rigid_Body*, Rigid_Bar_1*, Rigid_Bar_2*> bigArrayMaker(double *baseAngles, int numBasePendulums, int numPerturbations, double radius);










#endif /* HEADERS_DOUBLEPENDULUMCHAOS_H_ */
