
// Include header of Vector class because there are actual variables being defined
#include "vectors.h"

#ifndef MASS_H_
#define MASS_H_

#include "constraint_bodies.h"

// Class for a circular rigid body
class Circular_Rigid_Body
{
	public:

	// Constructors
	Circular_Rigid_Body();
	Circular_Rigid_Body(Vector &pos);
	Circular_Rigid_Body(double mass, Vector &pos, Vector &linear_vel, Vector &force_ext);

	// Position, Velocity
	Vector pos;
	Vector linear_vel;

	// External force
	Vector force_ext;

	// Mass
	double mass;

	// Find initial external force
	void set_initial_force_ext(void);

	// Copy data from to self from other Circular Rigid Body
	void copyRigidBody(Circular_Rigid_Body *o);
};



#endif /* MASS_H_ */
