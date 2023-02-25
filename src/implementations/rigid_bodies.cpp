#include "../headers/vectors.h"
#include "../headers/renderer.h"
#include "../headers/rigid_bodies.h"


// Class for a point mass
void Circular_Rigid_Body::set_initial_force_ext(void)
{
	// Initial external force is gravitational
	// Gravitational
	force_ext.y = mass * -9.8;
	force_ext.x = 0.0;
}

// Copy data from one Circular Rigid Body to another
void Circular_Rigid_Body::copyRigidBody(Circular_Rigid_Body *o)
{
	// Position and velocity
	pos = {o->pos.x, o->pos.y};
	linear_vel = {o->linear_vel.x, o->linear_vel.y};

	// External force
	force_ext = {o->force_ext.x, o->force_ext.y};

	// Mass
	mass = o->mass;
}





