
// Include header of Vector class because there are actual variables being defined
#include "vectors.h"

#ifndef MASS_H_
#define MASS_H_

#include "constraint_bodies.h"

// Class for a circular rigid body
class Circular_Rigid_Body
{
	public:

	// Position, Velocity
	Vector pos;
	Vector linear_vel;

	// External force
	Vector force_ext;

	// Mass
	double mass;

	// Find initial external force
	void set_initial_force_ext(void);
};



#endif /* MASS_H_ */
