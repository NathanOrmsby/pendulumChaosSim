/*
 * rk4.h
 *
 *  Created on: Nov 15, 2022
 *      Author: norms
 */

#ifndef RK4_H_
#define RK4_H_

#include "get_state.h"
#include "rigid_bodies.h"
#include "constraint_bodies.h"
#include "springs.h"

// Runge Kutta 4th Order to solve differential equations

// Optimized rk4
void rk4(Circular_Rigid_Body *mass_list, int numMasses, Rigid_Bar_1 *bar1s, int num_bar1s, Rigid_Bar_2 *bar2s, int num_bar2s, double dt);
void rk4_old(Circular_Rigid_Body *mass_list, int numMasses, Rigid_Bar_1 *bar1s, int num_bar1s, Rigid_Bar_2 *bar2s, int num_bar2s, double dt);

// Fills a mass list for another state
void copy_mass_list(Circular_Rigid_Body *old_list, Circular_Rigid_Body *new_list, int num_bodies);

#endif /* RK4_H_ */
