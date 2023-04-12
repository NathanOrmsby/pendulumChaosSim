// Calculates the kinetic plus potential energies of masses in the scene

#include <iostream>


#include "../headers/total_energy.h"

double total_energy(Circular_Rigid_Body *mass_list, int num_bodies)
{
	double tot_kinetic = 0.0;
	double tot_gp = 0.0;
	// Calculate the kinetic and potential energies of each mass

//	std::cout << "Loop count: " << loopCount << std::endl;
//	std::cout << "Num bodies: " << num_bodies << std::endl;
	for (int i = 0; i < num_bodies; i++)
	{
		double vx = mass_list[i].linear_vel.x;
		double vy = mass_list[i].linear_vel.y;
		tot_kinetic += 0.5 * mass_list[i].mass * ((vx * vx) + (vy * vy));
		tot_gp += mass_list[i].mass * 9.8 * mass_list[i].pos.y;

//		std::cout << "Body: " << i << std::endl;
//		std::cout << "Velocity: x: " << mass_list[i].linear_vel.x << " y: " << mass_list[i].linear_vel.y << std::endl;
//		std::cout << "Total kinetic: " << tot_kinetic << std::endl;
//		std::cout << "Total gravitational potential: " << tot_gp << std::endl;
	}
	return tot_kinetic + tot_gp;
}


