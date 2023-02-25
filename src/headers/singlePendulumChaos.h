/*
 * singlePendulumChaos.h
 *
 *  Created on: Feb 24, 2023
 *      Author: norms
 */

#ifndef HEADERS_SINGLEPENDULUMCHAOS_H_
#define HEADERS_SINGLEPENDULUMCHAOS_H_


// Simulates a single double pendulum with a circle of perturbations in the angle space
void singleDoublePendulumChaos(double *angles, double *barLens, double *massAmounts, double dt, int totalTime, int numPoints, double radius);


#endif /* HEADERS_SINGLEPENDULUMCHAOS_H_ */
