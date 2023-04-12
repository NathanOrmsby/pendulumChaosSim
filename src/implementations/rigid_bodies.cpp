#include "../headers/vectors.h"
#include "../headers/rigid_bodies.h"


// Constructors
Circular_Rigid_Body::Circular_Rigid_Body() : mass(0), pos({0, 0}), linear_vel({0, 0}), force_ext({0, 0}) {
}

Circular_Rigid_Body::Circular_Rigid_Body(Vector &pos) : mass(1.0), pos(pos), linear_vel({0, 0}), force_ext({0, 0}) {
}

Circular_Rigid_Body::Circular_Rigid_Body(double mass, Vector &pos, Vector &linear_vel, Vector &force_ext)
    : mass(mass), pos(pos), linear_vel(linear_vel), force_ext(force_ext) {
}

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





