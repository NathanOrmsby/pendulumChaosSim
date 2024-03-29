#include <cmath>
#include <iostream>

#include "../headers/constraint_bodies.h"

// The rigid bar classes. NEVER CHANGES LENGTH BETWEEN POINTS
// Rigid bar 1 methods
// This rigid bar attaches at a pivot at point0, and a mass at point1

// Constructors
Rigid_Bar_1::Rigid_Bar_1() : pivot({0, 0}), attached_mass(0) {
}

Rigid_Bar_1::Rigid_Bar_1(Vector &pivot, int attached_mass)
    : pivot(pivot), attached_mass(attached_mass) {
}

double Rigid_Bar_1::constraint(Circular_Rigid_Body *mass_list)
{
	// Given a state vector q, return a scalar value of this constraint function.
	// For single attachments, where one end is fixed. Connects a maximum of two masses
	// STATE VECTOR HAS 3 elements per mass, x, y, and theta

	// First constraint equation
	// Set the pivot to the origin
	// Shift x1 and y1 so they are relative to the pivot
	double x1 = mass_list[attached_mass].pos.x;
	double y1 = mass_list[attached_mass].pos.y;

	// Calculate the constraint
	double length = vector_magnitude({initial_point.x - pivot.x, initial_point.y - pivot.y});
	return ((x1 - pivot.x) * (x1 - pivot.x)) + ((y1 - pivot.y) * (y1 - pivot.y)) - (length * length);
}

double Rigid_Bar_1::constraint_time_derivative(Circular_Rigid_Body *mass_list)
{
	// Returns a value for the constraint time derivative
	double x1 = mass_list[attached_mass].pos.x;
	double y1 = mass_list[attached_mass].pos.y;
	double x1_dot = mass_list[attached_mass].linear_vel.x;
	double y1_dot = mass_list[attached_mass].linear_vel.y;
	return 2.0 * x1_dot * (x1 - pivot.x) + 2.0 * y1_dot * (y1 - pivot.y);

}

// Partial of constraint function with respect to x
double Rigid_Bar_1::jacobian_entry_x(Circular_Rigid_Body* mass_list)
{
	return 2.0 * (mass_list[attached_mass].pos.x - pivot.x);
}

// Partial of constraint function with respect to y
double Rigid_Bar_1::jacobian_entry_y(Circular_Rigid_Body *mass_list)
{
	return 2.0 * (mass_list[attached_mass].pos.y - pivot.y);
}
// partial of constraint time derivative wrt x and y
double Rigid_Bar_1::jacobian_derivative_entry_x(Circular_Rigid_Body *mass_list)
{
	return 2.0 * mass_list[attached_mass].linear_vel.x;
}

double Rigid_Bar_1::jacobian_derivative_entry_y(Circular_Rigid_Body *mass_list)
{
	return 2.0 * mass_list[attached_mass].linear_vel.y;
}

// Assign initial point given attached mass
void Rigid_Bar_1::determine_initial_point(Circular_Rigid_Body *mass_list)
{
	initial_point = mass_list[attached_mass].pos;
}

// Copy data from another Rigid_Bar_1 to self
void Rigid_Bar_1::copyRigidBar1(Rigid_Bar_1 *o)
{
	pivot = {o->pivot.x, o->pivot.y};
	initial_point = {o->initial_point.x, o->initial_point.y};
	attached_mass = o->attached_mass;
}

// Rigid bar 2 functions
// Attaches to two moving mass objects

// Constructor
Rigid_Bar_2::Rigid_Bar_2() : attached_masses{0, 1} {
}

Rigid_Bar_2::Rigid_Bar_2(int attached_mass_1, int attached_mass_2) {
    attached_masses[0] = attached_mass_1;
    attached_masses[1] = attached_mass_2;
}

double Rigid_Bar_2::constraint(Circular_Rigid_Body *mass_list)
{

		double x1 = mass_list[attached_masses[0]].pos.x;
		double y1 = mass_list[attached_masses[0]].pos.y;
		double x2 = mass_list[attached_masses[1]].pos.x;
		double y2 = mass_list[attached_masses[1]].pos.y;

		// Calculate the constraint
		double length = vector_magnitude({initial_point1.x - initial_point0.x, initial_point1.y - initial_point0.y});
		return ((x2 - x1) * (x2 - x1)) + ((y2 - y1) * (y2 - y1)) - (length * length);
}

double Rigid_Bar_2::constraint_time_derivative(Circular_Rigid_Body *mass_list)
{
	double x1 = mass_list[attached_masses[0]].pos.x;
	double y1 = mass_list[attached_masses[0]].pos.y;
	double x2 = mass_list[attached_masses[1]].pos.x;
	double y2 = mass_list[attached_masses[1]].pos.y;

	double x1_dot = mass_list[attached_masses[0]].linear_vel.x;
	double y1_dot = mass_list[attached_masses[0]].linear_vel.y;
	double x2_dot = mass_list[attached_masses[1]].linear_vel.x;
	double y2_dot = mass_list[attached_masses[1]].linear_vel.y;

	return 2.0 * x1_dot * (x1 - x2) + 2.0 * y1_dot * (y1 - y2) + 2.0 * x2_dot * (x2 - x1) + 2.0 * y2_dot * (y2 - y1);
}

double Rigid_Bar_2::jacobian_entry_x1(Circular_Rigid_Body *mass_list)
{

	double x1 = mass_list[attached_masses[0]].pos.x;
	double x2 = mass_list[attached_masses[1]].pos.x;
	return -2.0 * (x2 - x1);

}

double Rigid_Bar_2::jacobian_entry_y1(Circular_Rigid_Body *mass_list)
{

	double y1 = mass_list[attached_masses[0]].pos.y;
	double y2 = mass_list[attached_masses[1]].pos.y;
	return -2.0 * (y2 - y1);
}

double Rigid_Bar_2::jacobian_entry_x2(Circular_Rigid_Body *mass_list)
{

	double x1 = mass_list[attached_masses[0]].pos.x;
	double x2 = mass_list[attached_masses[1]].pos.x;
	return 2.0 * (x2 - x1);
}

double Rigid_Bar_2::jacobian_entry_y2(Circular_Rigid_Body *mass_list)
{

	double y1 = mass_list[attached_masses[0]].pos.y;
	double y2 = mass_list[attached_masses[1]].pos.y;
	return 2.0 * (y2 - y1);
}

// partial of constraint time derivative wrt x1, y2, x2, y2
double Rigid_Bar_2::jacobian_derivative_entry_x1(Circular_Rigid_Body *mass_list)
{
	double x1_dot = mass_list[attached_masses[0]].linear_vel.x;
	double x2_dot = mass_list[attached_masses[1]].linear_vel.x;
	return 2.0 * (x1_dot - x2_dot);
}

double Rigid_Bar_2::jacobian_derivative_entry_y1(Circular_Rigid_Body *mass_list)
{
	double y1_dot = mass_list[attached_masses[0]].linear_vel.y;
	double y2_dot = mass_list[attached_masses[1]].linear_vel.y;
	return 2.0 * (y1_dot - y2_dot);
}

double Rigid_Bar_2::jacobian_derivative_entry_x2(Circular_Rigid_Body *mass_list)
{
	double x1_dot = mass_list[attached_masses[0]].linear_vel.x;
	double x2_dot = mass_list[attached_masses[1]].linear_vel.x;
	return 2.0 * (x2_dot - x1_dot);
}

double Rigid_Bar_2::jacobian_derivative_entry_y2(Circular_Rigid_Body *mass_list)
{
	double y1_dot = mass_list[attached_masses[0]].linear_vel.y;
	double y2_dot = mass_list[attached_masses[1]].linear_vel.y;
	return 2.0 * (y2_dot - y1_dot);
}

void Rigid_Bar_2::determine_initial_points(Circular_Rigid_Body *mass_list)
{
	initial_point0 = mass_list[attached_masses[0]].pos;
	initial_point1 = mass_list[attached_masses[1]].pos;
}

// Copy data to self from other Rigid_Bar_2
void Rigid_Bar_2::copyRigidBar2(Rigid_Bar_2 *o)
{
	initial_point0 = {o->initial_point0.x, o->initial_point0.y};
	initial_point1 = {o->initial_point1.x, o->initial_point1.y};

	for (int i = 0; i < 2; ++i)
	{
		attached_masses[i] = o->attached_masses[i];
	}
}

