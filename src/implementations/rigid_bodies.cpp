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

	// Then apply spring force if a spring is attached using the spring class method
}





