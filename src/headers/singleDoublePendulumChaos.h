/*
 * singleDoublePendulumChaos.h
 *
 *  Created on: Apr 12, 2023
 *      Author: norms
 */

#ifndef HEADERS_SINGLEDOUBLEPENDULUMCHAOS_H_
#define HEADERS_SINGLEDOUBLEPENDULUMCHAOS_H_

// Simulation
double singleDoublePendulumChaosSim(double baseAngles[2], int numPerturbations, double perturbationRadius, double dt, int totalTime, int framesPerSecond, int lyapunovFrequency, bool data);

// Lyapunov calculation functions
// Timestep lyapunov function
void singleDoublePendulumTimestepLyapunov(Circular_Rigid_Body *masses, double *lyapunovRunningSums, int numPerturbations, double initialSeparation);

// Evaluate largest lyapunov exponent
double largestLyapunovResult(double *lyapunovRunningSums, int numPerturbations, const int totalSteps, const double dt);




#endif /* HEADERS_SINGLEDOUBLEPENDULUMCHAOS_H_ */
