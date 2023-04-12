/*
 * rk4.cpp
 *
 *  Created on: Nov 15, 2022
 *      Author: norms
 */


#include "../headers/rk4.h"


// Runge Kutta 4th Order to solve differential equations

// Optimized rk4
void rk4(Circular_Rigid_Body *mass_list, int numMasses, Rigid_Bar_1 *bar1s, int num_bar1s, Rigid_Bar_2 *bar2s, int num_bar2s, double dt)
{

	// Some preprocessing for optimization

	// Preprocess vector memory size
	int vectorLen = 2 * numMasses;
	int vectorMemSize = vectorLen * sizeof(double);

	// Malloc calls: Combine all forces and all velocities into single vector
	double *forces = (double *)malloc(5 * vectorMemSize);
	double *velocities = (double *)malloc(5 * vectorMemSize);

	// Memory indexing by k stage
	double *k1_force = forces;
	double *k1_vel = velocities;

	double *k2_force = forces + vectorLen;
	double *k2_vel = velocities + vectorLen;

	double *k3_force = forces + 2 * vectorLen;
	double *k3_vel = velocities + 2 * vectorLen;

	double *k4_force = forces + 3 * vectorLen;
	double *k4_vel = velocities + 3 * vectorLen;

	double *rk4_force = forces + 4 * vectorLen;
	double *rk4_vel = velocities + 4 * vectorLen;

	// k1 step
	// Create state: Use tmp state across k1-k4 to minimize unnecessary allocations
	State_Getter *tmp_state = new State_Getter(numMasses, num_bar1s + num_bar2s);

	// Get net force at t0
	// Calculate the net force at the starting position
	tmp_state->get_current_state(mass_list, bar1s, num_bar1s, bar2s, num_bar2s);

	// Copy net force list into k1_force, and populate k_vel
	for (int i = 0; i < numMasses; i++)
	{
	    // Index preprocessing
	    int ind1 = 2 * i;
	    int ind2 = ind1 + 1;

	    // Populate values
	    k1_vel[ind1] = mass_list[i].linear_vel.x;
	    k1_vel[ind2] = mass_list[i].linear_vel.y;
	    k1_force[ind1] = tmp_state->net_force_vector[ind1];
	    k1_force[ind2] = tmp_state->net_force_vector[ind2];
	}

	// Free the internal arrays within State_Getter
	tmp_state->~State_Getter();

	// k2 step
	// New mass_list for k2 state: Use a tmp mass list that can be repurposed.
	Circular_Rigid_Body *tmp_mass_list = (Circular_Rigid_Body *)malloc(numMasses * sizeof(Circular_Rigid_Body));

	// Copy properties of mass_list
	copy_mass_list(mass_list, tmp_mass_list, numMasses);


	// Timestep preprocessing
	double half_step = 0.5 * dt;
	// Push k2 state forward dt / 2 using k1 as net force
	for (int i = 0; i < numMasses; i++)
	{
		// Index preprocessing
		int ind1 = 2 * i;
		int ind2 = ind1 + 1;

		// Pointer preprocessing
		Circular_Rigid_Body &mass = tmp_mass_list[i];

		// Value preprocessing
		double accelx = k1_force[ind1] / mass.mass;
		double accely = k1_force[ind2] / mass.mass;

		// Update velocities:
		mass.linear_vel.x += accelx * half_step;
		mass.linear_vel.y += accely * half_step;

		// PREVIOUS
		// Update positions: Using k1 velocity. Initial position velocity
		mass.pos.x += k1_vel[ind1] * half_step;
		mass.pos.y += k1_vel[ind2] * half_step;
	}

	// State_getter for k2: Reuse tmp
	tmp_state->get_current_state(tmp_mass_list, bar1s, num_bar1s, bar2s, num_bar2s);

	// Copy net force list into k2_force, and populate k_vel
	for (int i = 0; i < numMasses; i++)
	{
		// Index preprocessing
		int ind1 = 2 * i;
		int ind2 = ind1 + 1;

		// Populate values
		k2_vel[ind1] = tmp_mass_list[i].linear_vel.x;
		k2_vel[ind2] = tmp_mass_list[i].linear_vel.y;
		k2_force[ind1] = tmp_state->net_force_vector[ind1];
		k2_force[ind2] = tmp_state->net_force_vector[ind2];
	}


	// Free the internal arrays within State_Getter: BROKEN RIGHT HERE
	tmp_state->~State_Getter();

	// k3 step
	// New mass_list for k3: use tmp
	copy_mass_list(mass_list, tmp_mass_list, numMasses);

	// Push k3 state forward dt / 2 using k1 as net force
	for (int i = 0; i < numMasses; i++)
	{
		// Index preprocessing
		int ind1 = 2 * i;
		int ind2 = ind1 + 1;

		// Pointer preprocessing
		Circular_Rigid_Body &mass = tmp_mass_list[i];

		// Value preprocessing
		double accelx = k2_force[ind1] / mass.mass;
		double accely = k2_force[ind2] / mass.mass;

		// Update velocities:
		mass.linear_vel.x += accelx * half_step;
		mass.linear_vel.y += accely * half_step;

		// PREVIOUS
		// Update positions: Using k1 velocity. Initial position velocity
		mass.pos.x += k2_vel[ind1] * half_step;
		mass.pos.y += k2_vel[ind2] * half_step;
	}

	// New state_getter for k3

	// Get state for k3
	tmp_state->get_current_state(tmp_mass_list, bar1s, num_bar1s, bar2s, num_bar2s);

	// Copy net force list into k3_force, and populate k3_vel
	for (int i = 0; i < numMasses; i++)
	{
		// Index preprocessing
		int ind1 = 2 * i;
		int ind2 = ind1 + 1;

		// Populate values
		k3_vel[ind1] = tmp_mass_list[i].linear_vel.x;
		k3_vel[ind2] = tmp_mass_list[i].linear_vel.y;
		k3_force[ind1] = tmp_state->net_force_vector[ind1];
		k3_force[ind2] = tmp_state->net_force_vector[ind2];
	}

	// Free the internal arrays within State_Getter
	tmp_state->~State_Getter();

	// k4 step
	// Make new mass_list for k4 state
	copy_mass_list(mass_list, tmp_mass_list, numMasses);

	// Push k4 state forward dt using k3 force and vel
	for (int i = 0; i < numMasses; i++)
	{
		// Index preprocessing
		int ind1 = 2 * i;
		int ind2 = ind1 + 1;

		// Pointer preprocessing
		Circular_Rigid_Body &mass = tmp_mass_list[i];

		// Value preprocessing
		double accelx = k3_force[ind1] / mass.mass;
		double accely = k3_force[ind2] / mass.mass;

		// Update velocities:
		mass.linear_vel.x += accelx * dt;
		mass.linear_vel.y += accely * dt;

		// PREVIOUS
		// Update positions: Using k1 velocity. Initial position velocity
		mass.pos.x += k3_vel[ind1] * dt;
		mass.pos.y += k3_vel[ind2] * dt;
	}

	// New state_getter for k4

	// Get state for k4
	tmp_state->get_current_state(tmp_mass_list, bar1s, num_bar1s, bar2s, num_bar2s);

	// Copy net force list into k4_force, and populate k4_vel
	for (int i = 0; i < numMasses; i++)
	{
		// Index preprocessing
		int ind1 = 2 * i;
		int ind2 = ind1 + 1;

		// Populate values
		k4_vel[ind1] = tmp_mass_list[i].linear_vel.x;
		k4_vel[ind2] = tmp_mass_list[i].linear_vel.y;
		k4_force[ind1] = tmp_state->net_force_vector[ind1];
		k4_force[ind2] = tmp_state->net_force_vector[ind2];
	}

	// Calculate the weighted rk4 force and velocity
	for (int i = 0; i < vectorLen; i++)
	{
		rk4_force[i] = k1_force[i] + 2 * k2_force[i] + 2 * k3_force[i] + k4_force[i];
		rk4_vel[i] = k1_vel[i] + 2 * k2_vel[i] + 2 * k3_vel[i] + k4_vel[i];
	}

	// Push the initial state forward dt / 6 using the weighted average rk4 force as net force

	// Timestep preprocessing
	double one_sixth_step = dt / 6.0;

	for (int i = 0; i < numMasses; i++)
	{
		// Index preprocessing
		int ind1 = 2 * i;
		int ind2 = ind1 + 1;

		// Pointer preprocessing
		Circular_Rigid_Body &mass = mass_list[i];

		// Value preprocessing
		double accelx = rk4_force[ind1] / mass.mass;
		double accely = rk4_force[ind2] / mass.mass;

		// Update velocities:
		mass.linear_vel.x += accelx * one_sixth_step;
		mass.linear_vel.y += accely * one_sixth_step;

		// PREVIOUS
		// Update positions: Using k1 velocity. Initial position velocity
		mass.pos.x += rk4_vel[ind1] * one_sixth_step;
		mass.pos.y += rk4_vel[ind2] * one_sixth_step;
	}

	// Free all the stuff
	free(forces);
	free(velocities);
	free(tmp_mass_list);

	// Delete tmp_state, calls destructor
	delete tmp_state;

}

// Pre optimization version: Made a few changes, if it doesn't work just swipe from another version that worked.
void rk4_old(Circular_Rigid_Body *mass_list, int numMasses, Rigid_Bar_1 *bar1s, int num_bar1s, Rigid_Bar_2 *bar2s, int num_bar2s, double dt)
{
	// Some preprocessing for optimization

	// Preprocess vector memory size
	int vectorLen = 2 * numMasses;
	int vectorMemSize = vectorLen * sizeof(double);

	// k1 step
	// Create state
	State_Getter *k1_state = new State_Getter(numMasses, num_bar1s + num_bar2s);

	// Get net force at t0
	// Calculate the net force at the starting position
	k1_state->get_current_state(mass_list, bar1s, num_bar1s, bar2s, num_bar2s);

	// Copy net force list into k1
	double *k1_force = (double *)malloc(vectorMemSize);
	for (int i = 0; i < vectorLen; i++)
	{
		k1_force[i] = k1_state->net_force_vector[i];
	}

	delete k1_state;


	// K1 velocity
	double *k1_vel = (double *)malloc(vectorMemSize);
	for (int i = 0; i < numMasses; i++)
	{
		k1_vel[2 * i] = mass_list[i].linear_vel.x;
		k1_vel[2 * i + 1] = mass_list[i].linear_vel.y;
	}

	// k2 step
	// Make new mass_list for k2 state.
	Circular_Rigid_Body *tmp_mass_list = (Circular_Rigid_Body *)malloc(numMasses * sizeof(Circular_Rigid_Body));

	// Copy properties of mass_list
	copy_mass_list(mass_list, tmp_mass_list, numMasses);

	// Push k2 state forward dt / 2 using k1 as net force
	for (int i = 0; i < numMasses; i++)
	{
		// Update velocities: THIS USES k1 / 2
		tmp_mass_list[i].linear_vel.x += (k1_force[2 * i] / tmp_mass_list[i].mass) * 0.5 * dt;
		tmp_mass_list[i].linear_vel.y += (k1_force[2 * i + 1] / tmp_mass_list[i].mass * 0.5) * dt;

		// PREVIOUS
		// Update positions: Using k1 velocity. Initial position velocity
		tmp_mass_list[i].pos.x += k1_vel[2 * i] * 0.5 * dt;
		tmp_mass_list[i].pos.y += k1_vel[2 * i + 1] * 0.5 * dt;
	}

	// New state_getter for k2
	State_Getter *k2_state = (State_Getter *)malloc(sizeof(State_Getter));
	k2_state->num_bodies = numMasses;
	k2_state->num_constraints = num_bar1s + num_bar2s;

	// Get state for k2
	k2_state->get_current_state(tmp_mass_list, bar1s, num_bar1s, bar2s, num_bar2s);

	// Fill up k2 net force
	double *k2_force = (double *)malloc(vectorMemSize);
	for (int i = 0; i < vectorLen; i++)
	{
		k2_force[i] = k2_state->net_force_vector[i];
	}

	// Get k2 velocity vector
	double *k2_vel = (double *)malloc(vectorMemSize);
	for (int i = 0; i < numMasses; i++)
	{

		k2_vel[2 * i] = tmp_mass_list[i].linear_vel.x;
		k2_vel[2 * i + 1] = tmp_mass_list[i].linear_vel.y;
	}

	// Free k2 mass list
	free(tmp_mass_list);
	// Free k2 state
	delete k2_state;

	// k3 step
	// Make new mass_list for k3 state
	Circular_Rigid_Body *mass_list_k3 = (Circular_Rigid_Body *)malloc(numMasses * sizeof(Circular_Rigid_Body));
	copy_mass_list(mass_list, mass_list_k3, numMasses);

	// Push the mass_list_k3 state forward dt / 2 using k2 as net force
	for (int i = 0; i < numMasses; i++)
	{
		// Update velocities
		mass_list_k3[i].linear_vel.x += (k2_force[2 * i] / mass_list_k3[i].mass) * 0.5 * dt;
		mass_list_k3[i].linear_vel.y += (k2_force[2 * i + 1] / mass_list_k3[i].mass) * 0.5 * dt;

		// PREVIOUS
		// Update positions
		mass_list_k3[i].pos.x += k2_vel[2 * i] * 0.5 * dt;
		mass_list_k3[i].pos.y += k2_vel[2 * i + 1] * 0.5 * dt;
	}

	// New state_getter for k3
	State_Getter *k3_state = (State_Getter *)malloc(sizeof(State_Getter));
	k3_state->num_bodies = numMasses;
	k3_state->num_constraints = num_bar1s + num_bar2s;

	// Get state for k3
	k3_state->get_current_state(mass_list_k3, bar1s, num_bar1s, bar2s, num_bar2s);

	// Fill up k3 net force
	double *k3_force = (double *)malloc(vectorMemSize);
	for (int i = 0; i < vectorLen; i++)
	{
		k3_force[i] = k3_state->net_force_vector[i];
	}

	// Get k3 velocity vector
	double *k3_vel = (double *)malloc(vectorMemSize);
	for (int i = 0; i < numMasses; i++)
	{
		k3_vel[2 * i] = mass_list_k3[i].linear_vel.x;
		k3_vel[2 * i + 1] = mass_list_k3[i].linear_vel.y;
	}

	// Free stuff
	free(mass_list_k3);
	delete k3_state;

	// k4 step
	// Make new mass_list for k4 state

	Circular_Rigid_Body *mass_list_k4 = (Circular_Rigid_Body *)malloc(numMasses * sizeof(Circular_Rigid_Body));
	copy_mass_list(mass_list, mass_list_k4, numMasses);

	// Push the k4 state forward dt using k3 as net force and k3 vel
	for (int i = 0; i < numMasses; i++)
	{
		// Update velocities
		mass_list_k4[i].linear_vel.x += (k3_force[2 * i] / mass_list_k4[i].mass) * dt;
		mass_list_k4[i].linear_vel.y += (k3_force[2 * i + 1] / mass_list_k4[i].mass) * dt;

		// PREVIOUS
		// Update positions
		mass_list_k4[i].pos.x += k3_vel[2 * i] * dt;
		mass_list_k4[i].pos.y += k3_vel[2 * i + 1] * dt;
	}

	// New state_getter for k4
	State_Getter *k4_state = (State_Getter *)malloc(sizeof(State_Getter));
	k4_state->num_bodies = numMasses;
	k4_state->num_constraints = num_bar1s + num_bar2s;

	// Get state for k4
	k4_state->get_current_state(mass_list_k4, bar1s, num_bar1s, bar2s, num_bar2s);

	// Fill up k4 net force
	double *k4_force = (double *)malloc(vectorMemSize);
	for (int i = 0; i < vectorLen; i++)
	{
		k4_force[i] = k4_state->net_force_vector[i];
	}

	// Get k4 velocity vector
	double *k4_vel = (double *)malloc(vectorMemSize);
	for (int i = 0; i < numMasses; i++)
	{
		k4_vel[2 * i] = mass_list_k4[i].linear_vel.x;
		k4_vel[2 * i + 1] = mass_list_k4[i].linear_vel.y;
	}

	// Free stuff
	free(mass_list_k4);
	delete k4_state;

	// Calculate the weighted rk4 net force (k1_force + 2*k2_force + 2*k3_force + k4_force)
	double *rk4_force = (double *)malloc(vectorMemSize);
	double *rk4_vel = (double *)malloc(vectorMemSize);
	for (int i = 0; i < vectorLen; i++)
	{
		rk4_force[i] = k1_force[i] + 2 * k2_force[i] + 2 * k3_force[i] + k4_force[i];
	}
	for (int i = 0; i < vectorLen; i++)
	{
		rk4_vel[i] = k1_vel[i] + 2 * k2_vel[i] + 2 * k3_vel[i] + k4_vel[i];
	}

	// Calculate the weighted rk4 velocity (k1_vel + 2 * k2_vel + 2 * k3_vel + k4_vel)

	// Push the initial state forward dt using the weighted average rk4 force as net force
	for (int i = 0; i < numMasses; i++)
	{
		// Update velocities
		mass_list[i].linear_vel.x += ((rk4_force[2 * i] / mass_list[i].mass) / 6.0) * dt;
		mass_list[i].linear_vel.y += ((rk4_force[2 * i + 1] / mass_list[i].mass) / 6.0) * dt;

		// Update positions
		mass_list[i].pos.x += (rk4_vel[2 * i] / 6.0)  * dt;
		mass_list[i].pos.y += (rk4_vel[2 * i + 1] / 6.0)  * dt;
	}

	// Free all the stuff
	free(k1_force);
	free(k1_vel);
	free(k2_force);
	free(k2_vel);
	free(k3_force);
	free(k3_vel);
	free(k4_force);
	free(k4_vel);
	free(rk4_force);
	free(rk4_vel);
}

// Fills a mass list for another state
void copy_mass_list(Circular_Rigid_Body *old_list, Circular_Rigid_Body *new_list, int num_bodies)
{
	for (int i = 0; i < num_bodies; i++)
	{
		new_list[i].pos.x = old_list[i].pos.x;
		new_list[i].pos.y = old_list[i].pos.y;
		new_list[i].linear_vel.x = old_list[i].linear_vel.x;
		new_list[i].linear_vel.y = old_list[i].linear_vel.y;
		new_list[i].force_ext.x = old_list[i].force_ext.x;
		new_list[i].force_ext.y = old_list[i].force_ext.y;
		new_list[i].mass = old_list[i].mass;
	}
}
