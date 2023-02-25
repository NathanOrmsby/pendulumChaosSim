/*
 * Pendulum.cpp
 *
 *  Created on: Feb 24, 2023
 *      Author: norms
 */

#include "../headers/Pendulum.h"
#include "Vectors.h"

#include <math.h>

// Methods that act on the Pendulum Struct

// Initialize a Pendulum
Pendulum *initializePendulum(Circular_Rigid_Body *masses, int numMasses, Rigid_Bar_1 bar1, Rigid_Bar_2 *bar2s, int numBar2s)
{
	// Allocate space
	Pendulum *p = (Pendulum *)malloc(sizeof(Pendulum));

	// Set the masses
	p->masses = masses;
	p->numMasses = numMasses;

	// Set the pivot bar
	p->bar1 = bar1;

	// Set the bars connecting the masses if there are any
	p->numBar2s = numBar2s;
	if (p->numBar2s > 0)
	{
		p->bar2s = bar2s;
	}

	// Set the angles for the Pendulum
	p->numAngles = numMasses;
	p->angles = (double *)malloc(p->numAngles * sizeof(double));

	// Set pendulum angles
	p->angles = findPendulumAngles(p);

	return p;
}

// Initialize a double pendulum centered at the origin
Pendulum *initializeDoublePendulum(double *angles, double *barLens, double *massAmounts)
{
	// Allocate the pendulum and all fields
	Pendulum *p = (Pendulum *)malloc(sizeof(Pendulum));
	p->masses = (Circular_Rigid_Body *)malloc(sizeof(Circular_Rigid_Body));
	p->angles = (double *)malloc(2 * sizeof(double));
	p->bar2s = (Rigid_Bar_2 *)malloc(sizeof(Rigid_Bar_2));

	// Trivial initializations
	p->numAngles = 2;
	p->numMasses = 2;
	p->numBar2s = 1;
	p->perturbed = nullptr;

	// Angles
	for (int i = 0; i < 2; ++i)
	{
		p->angles[i] = angles[i];
	}

	// Initialize positions and mass amounts of masses
	p->masses = initializeMassPositions(angles, p->numAngles, barLens, p->numBar2s + 1, massAmounts);

	// Further mass initializations
	for (int i = 0; i < p->numMasses; ++i)
	{
		// Set velocity to zero
		p->masses[i].linear_vel = {0, 0};

		// Set initial gravitational force
		p->masses[i].set_initial_force_ext();
	}

	// Initialize the rigid bars

	// Bar1
	p->bar1.attached_mass = 0;
	p->bar1.pivot = {0, 0};
	p->bar1.determine_initial_point(p->masses);

	// Bar2
	p->bar2s[0].attached_masses = {0, 1};
	p->bar2s[0].determine_initial_points(p->masses);

	return p;
}

// Initialize mass positions given angles and bar lengths assuming centered at origin
Circular_Rigid_Body *initializeMassPositions(double *angles, int numAngles, double *barLens, int numBars, double *massAmounts)
{
	// Number of masses is equivalent to number of angles
	int numMasses = numAngles;

	// Allocate space
	Circular_Rigid_Body *masses = (Circular_Rigid_Body *)malloc(numMasses * sizeof(Circular_Rigid_Body));

	// Find cartesian coordinates for all masses
	// Running offset for shifted origin
	Vector offset = {0, 0};

	for (int i = 0; i < numMasses; ++i)
	{
		// Set mass
		masses[i].mass = massAmounts[i];

		// Calculate position
		masses[i].pos.x = barLens[i] * sin(angles[i]) + offset.x;
		masses[i].pos.y = barLens[i] * cos(angles[i]) + offset.y;

		// Update the offset
		offset.x += masses[i].pos.x;
		offset.y += masses[i].pos.y;
	}

	return masses;
}

// Return a double array of all angles for a pendulum
double *findPendulumAngles(Pendulum *p)
{
	// Allocate for the result
	double *angles = (double *)malloc(p->numAngles * sizeof(double));

	// Angle for bar between pivot and first mass
	Vector v = {p->masses[0].pos.x - p->bar1.pivot.x, p->masses[0].pos.y - p->bar1.pivot.y};
	Vector vert = {0, -1};
	angles[0] = angle_between_vectors(v, vert);

	// Find the rest of the angles
	if (p->numAngles > 1)
	{
		for (int i = 1; i < p->numMasses; ++i)
		{
			v = {p->masses[i].pos.x - p->masses[i - 1].pos.x, p->masses[i].pos.y - p->masses[i - 1].pos.y};
			angles[i] = angle_between_vectors(v, vert);
		}
	}

	return angles;
}

// Copies all pendulum data from i to o
void copyPendulum(Pendulum *o, Pendulum *i)
{
	o->numMasses = i->numMasses;
	o->numBar2s = i->numBar2s;
	o->numAngles = i->numAngles;

	// Copy mass objects

	// Allocate space
	o->masses = (Circular_Rigid_Body *)malloc(o->numMasses * sizeof(Circular_Rigid_Body));

	// Copy data
	for (int j = 0; j < o->numMasses; ++j)
	{
		o->masses[j].copyRigidBody(&i->masses[j]);
	}

	// Copy rigid bars
	o->bar1.copyRigidBar1(&i->bar1);

	// Allocate space for bar2s
	o->bar2s = (Rigid_Bar_2 *)malloc(o->numBar2s * sizeof(Rigid_Bar_2));

	// Copy data
	for (int j = 0; j < o->numBar2s; ++j)
	{
		o->bar2s[j].copyRigidBar2(&i->bar2s[j]);
	}
}

// Free a pendulum struct from memory
void freePendulum(Pendulum *p)
{
	free(p->masses);
	free(p->bar2s);
	free(p->angles);

	if (p->perturbed != nullptr)
	{
		free(p->perturbed);
	}

	free(p);
}




