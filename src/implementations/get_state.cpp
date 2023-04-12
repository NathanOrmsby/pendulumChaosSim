#include <iostream>
#include <chrono>


#include "../headers/get_state.h"
#include "../headers/bc_gradient.h"
#include "../headers/matrix_stuff.h"

// New constructors for new matrix block style
// Default
State_Getter::State_Getter() : num_bodies(0), num_constraints(0), state_vector(nullptr), state_vector_derivative(nullptr),
                        inverse_mass_matrix(nullptr), force_ext_vector(nullptr), constraint_vector(nullptr),
                        constraint_derivative_vector(nullptr), jacobian(nullptr), jacobian_data(nullptr), jacobian_derivative(nullptr),
						jacobian_derivative_data(nullptr), constraint_force_vector(nullptr), net_force_vector(nullptr),ks(0.0), kd(0.0) {}

// Specifying non pointer attributes
State_Getter::State_Getter(int num_bodies, int num_constraints) :
    num_bodies(num_bodies),
    num_constraints(num_constraints),
    state_vector(nullptr),
    state_vector_derivative(nullptr),
    inverse_mass_matrix(nullptr),
    force_ext_vector(nullptr),
    constraint_vector(nullptr),
    constraint_derivative_vector(nullptr),
    jacobian(nullptr),
	jacobian_data(nullptr),
    jacobian_derivative(nullptr),
	jacobian_derivative_data(nullptr),
    constraint_force_vector(nullptr),
    net_force_vector(nullptr) {}

// Destructor
State_Getter::~State_Getter() {
	// Free the vectors
	free(state_vector);
	free(state_vector_derivative);
	free(inverse_mass_matrix);
	free(force_ext_vector);
	free(constraint_vector);
	free(constraint_derivative_vector);
	free(constraint_force_vector);
	free(net_force_vector);
	free(jacobian);
	free(jacobian_data);
	free(jacobian_derivative);
	free(jacobian_derivative_data);
}

//// OLD CONSTRUCTORS and Destructor FOR WORKING MATRIX BLOCK STATE_GETTER CLASS
//// Default
//State_Getter::State_Getter() : num_bodies(0), num_constraints(0), state_vector(nullptr), state_vector_derivative(nullptr),
//                        inverse_mass_matrix(nullptr), force_ext_vector(nullptr), constraint_vector(nullptr),
//                        constraint_derivative_vector(nullptr), jacobian(nullptr), jacobian_derivative(nullptr), constraint_force_vector(nullptr),
//						net_force_vector(nullptr),ks(0.0), kd(0.0) {}
//
//// Specifying non pointer attributes
//State_Getter::State_Getter(int num_bodies, int num_constraints) :
//    num_bodies(num_bodies),
//    num_constraints(num_constraints),
//    state_vector(nullptr),
//    state_vector_derivative(nullptr),
//    inverse_mass_matrix(nullptr),
//    force_ext_vector(nullptr),
//    constraint_vector(nullptr),
//    constraint_derivative_vector(nullptr),
//    jacobian(nullptr),
//    jacobian_derivative(nullptr),
//    constraint_force_vector(nullptr),
//    net_force_vector(nullptr) {}
//
//// Destructor
//State_Getter::~State_Getter() {
//	// Free the vectors
//	free(state_vector);
//	free(state_vector_derivative);
//	free(inverse_mass_matrix);
//	free(force_ext_vector);
//	free(constraint_vector);
//	free(constraint_derivative_vector);
//	free(constraint_force_vector);
//	free(net_force_vector);
//
//	// Free the matrix block structs
//	for (int i = 0; i < num_constraints; i++) {
//	// Free the matrix arrays in each block
//	free(jacobian[i].matrix);
//	}
//	free(jacobian);
//
//	for (int i = 0; i < num_constraints; i++) {
//	// Free the matrix arrays in each block
//	free(jacobian_derivative[i].matrix);
//	}
//	free(jacobian_derivative);
//}
// Functions involved with calculating the state of the system


// New version for new Matrix Block Class
void State_Getter::get_current_state(Circular_Rigid_Body *mass_list, Rigid_Bar_1 *bar1s, int num_bar1s, Rigid_Bar_2 *bar2s, int num_bar2s)
{
	// Results in a net force vector for the current state assigned to the class.
	// Net force vector is the sum of the external force vector and the constraint force vector.

	// *** Creating the matrices to start ***
	// Create state vector and state vector derivative

	create_state_vector(mass_list);
	create_state_vector_derivative(mass_list);
	// Create inverse mass matrix
	create_inverse_mass_matrix(mass_list);
	// Create external force vector
	create_force_ext_vector(mass_list);
	// Create jacobian
	create_jacobian(mass_list, bar1s, num_bar1s, bar2s, num_bar2s);
	// Create jacobian derivative
	create_jacobian_derivative(mass_list, bar1s, num_bar1s, bar2s, num_bar2s);
	// For spring and damping coefficients to remove numerical drift
	// Create constraint vector
	create_constraint_vector(mass_list, bar1s, num_bar1s, bar2s, num_bar2s);
	// Create constraint derivative vector
	create_constraint_derivative_vector(mass_list, bar1s, num_bar1s, bar2s, num_bar2s);



	// Now calculate b.
	double *b = new double[num_constraints]();
	calculate_b(b);


	// Solve for the lambda vector using biconjugate gradient method. try itol 1 to start
	double *x = new double[num_constraints]();

	// Try putting the max number of iterations at 10 * N
	// Iteration count and error in result.
	int *iter = (int *)malloc(sizeof(int));
	double *err = (double *)malloc(sizeof(double));

	biconjugate_gradient(this, b, x, 1, 1e-12, 10 * num_constraints, iter, err);


	// Free dynamically allocated iter and err.
	free(iter);
	free(err);

	// Free the rhs b
	delete[] b;


	// Now we have x, which is our lambda_vector.
	// Multiply by jacobian transpose to get constraint forces. F = J_T * lambda
	// Allocate space
	constraint_force_vector = new double[2 * num_bodies]();

	matrix_blocks_transpose_times_vector(jacobian, num_constraints, jacobian_data, x, constraint_force_vector);

	// Find the net force vector
	calculate_net_force_vector();

	// Free the minimized x vector
	delete[] x;
}

// WORKING VERSION WITH OLD MATRIX BLOCK CLASS
//void State_Getter::get_current_state(Circular_Rigid_Body *mass_list, Rigid_Bar_1 *bar1s, int num_bar1s, Rigid_Bar_2 *bar2s, int num_bar2s)
//{
//	// Results in a net force vector for the current state assigned to the class.
//	// Net force vector is the sum of the external force vector and the constraint force vector.
//
//	// *** Creating the matrices to start ***
//	// Create state vector and state vector derivative
//
//	create_state_vector(mass_list);
//	create_state_vector_derivative(mass_list);
//	// Create inverse mass matrix
//	create_inverse_mass_matrix(mass_list);
//	// Create external force vector
//	create_force_ext_vector(mass_list);
//	// Create jacobian
//	create_jacobian(mass_list, bar1s, num_bar1s, bar2s, num_bar2s);
//	// Create jacobian derivative
//	create_jacobian_derivative(mass_list, bar1s, num_bar1s, bar2s, num_bar2s);
//	// For spring and damping coefficients to remove numerical drift
//	// Create constraint vector
//	create_constraint_vector(mass_list, bar1s, num_bar1s, bar2s, num_bar2s);
//	// Create constraint derivative vector
//	create_constraint_derivative_vector(mass_list, bar1s, num_bar1s, bar2s, num_bar2s);
//
//
//
//	// Now calculate b.
//	double *b = new double[num_constraints]();
//	calculate_b(b);
//
//
//	// Solve for the lambda vector using biconjugate gradient method. try itol 1 to start
//	double *x = new double[num_constraints]();
//
//	// Try putting the max number of iterations at 10 * N
//	// Iteration count and error in result.
//	int *iter = (int *)malloc(sizeof(int));
//	double *err = (double *)malloc(sizeof(double));
//
//	biconjugate_gradient(this, b, x, 1, 1e-12, 10 * num_constraints, iter, err);
//
//
//	// Free dynamically allocated iter and err.
//	free(iter);
//	free(err);
//
//	// Free the rhs b
//	delete[] b;
//
//
//	// Now we have x, which is our lambda_vector.
//	// Multiply by jacobian transpose to get constraint forces. F = J_T * lambda
//	// Allocate space
//	constraint_force_vector = new double[2 * num_bodies]();
//
//	matrix_blocks_transpose_times_vector(jacobian, num_constraints, x, constraint_force_vector);
//
//	// Find the net force vector
//	calculate_net_force_vector();
//
//	// Free the minimized x vector
//	delete[] x;
//}

// Create the state_vector
void State_Getter::create_state_vector(Circular_Rigid_Body *mass_list)
{
	// Dynamically allocate memory
	state_vector = (double *)malloc(2 * num_bodies * sizeof(double));
	// Loop through the rigid bodies
	for (int i = 0; i < num_bodies; i++)
	{
		// x position
		state_vector[2 * i] = mass_list[i].pos.x;
		// y position
		state_vector[2 * i + 1] = mass_list[i].pos.y;
	}
}

// Creating the matrices and vectors
void State_Getter::create_state_vector_derivative(Circular_Rigid_Body *mass_list)
{
	// Dynamically allocate memory
	state_vector_derivative = (double *)malloc(2 * num_bodies * sizeof(double));
	// Loop through the rigid bodies
	for (int i = 0; i < num_bodies; i++)
	{
		// x position
		state_vector_derivative[2 * i] = mass_list[i].linear_vel.x;
		// y position
		state_vector_derivative[2 * i + 1] = mass_list[i].linear_vel.y;
	}
}

void State_Getter::create_inverse_mass_matrix(Circular_Rigid_Body *mass_list)
{
	// Dynamically allocate memory
	inverse_mass_matrix = (double *)malloc(2 * num_bodies * sizeof(double));
	// Loop through the rigid bodies
	for (int i = 0; i < num_bodies; i++)
	{
		// x position
		inverse_mass_matrix[2 * i] = 1.0 / mass_list[i].mass;
		// y position
		inverse_mass_matrix[2 * i + 1] = 1.0 / mass_list[i].mass;
	}
}

void State_Getter::create_force_ext_vector(Circular_Rigid_Body *mass_list)
{
	// Allocate space for the external force vector
	force_ext_vector = (double *)malloc(2 * num_bodies * sizeof(double));
	// Loop through masses and assign force vectors
	for (int i = 0; i < num_bodies; i++)
	{
		// Determine the external force on each mass
		mass_list[i].set_initial_force_ext();
		// Yoink the numbers from each mass: fx, fy, torque
		force_ext_vector[2 * i] = mass_list[i].force_ext.x;
		force_ext_vector[2 * i + 1] = mass_list[i].force_ext.y;
	}
}

void State_Getter::create_constraint_vector(Circular_Rigid_Body *mass_list, Rigid_Bar_1 *bar1s, int num_bar1s, Rigid_Bar_2 *bar2s, int num_bar2s)
{
	// Vector of scalar constraint values for each constraint
	// Allocate space
	constraint_vector = (double *)malloc(num_constraints * sizeof(double));
	// Loop through lists of each type of constraint, these are all the constraint functions
	// Rigid_Bar_1s.
	for (int i = 0; i < num_bar1s; i++)
	{
		// Each rigid_bar_1 applies a constraint to the vector
		constraint_vector[i] = bar1s[i].constraint(mass_list);
	}
	// Rigid_Bar_2s
	if (num_bar2s > 0)
	{
		for (int i = 0; i < num_bar2s; i++)
		{
			// Each rigid_bar_2 applies a constraint to the vector
			constraint_vector[num_bar1s + i] = bar2s[i].constraint(mass_list);
		}
	}
}

void State_Getter::create_constraint_derivative_vector(Circular_Rigid_Body *mass_list, Rigid_Bar_1 *bar1s, int num_bar1s, Rigid_Bar_2 *bar2s, int num_bar2s)
{
	// Vector of constraint time derivative values for each constraint
	// Allocate space
	constraint_derivative_vector = (double *)malloc(num_constraints * sizeof(double));

	// Loop through lists of each type of constraint, these are all the constraint functions
	// Rigid_Bar_1s.
	for (int i = 0; i < num_bar1s; i++)
	{
		// Each rigid_bar_1 applies a constraint to the vector
		constraint_derivative_vector[i] = bar1s[i].constraint_time_derivative(mass_list);
	}
	if (num_bar2s > 0)
	{
		// Rigid_Bar_2s
		for (int i = 0; i < num_bar2s; i++)
		{
			// Each rigid_bar_2 applies a constraint to the vector
			constraint_derivative_vector[num_bar1s + i] = bar2s[i].constraint_time_derivative(mass_list);
		}
	}
}

// Testing new create jacobian with new matrix block style: NEEDS TESTING
void State_Getter::create_jacobian(Circular_Rigid_Body *mass_list, Rigid_Bar_1 *bar1s, int num_bar1s, Rigid_Bar_2 *bar2s, int num_bar2s)
{
    // Allocate memory for the matrix block list, one block for each constraint
    jacobian = (Matrix_Block *)malloc(num_constraints * sizeof(Matrix_Block));

    // Allocate memory for the entire Jacobian array at once
    jacobian_data = (double *)malloc((2 * num_bar1s + num_bar2s * 4) * sizeof(double));

    // Loop through each rigid_bar_1 constraint
    for (int i = 0; i < num_bar1s; i++)
    {
        // Calculate the index in the flattened Jacobian array for the current constraint
        int jacobian_index = 2 * i;

        // Create a Matrix_Block for the current constraint and store its information
        jacobian[i] = Matrix_Block(1, 2, i, 2 * bar1s[i].attached_mass, jacobian_index);

        // Compute the Jacobian entries for the current constraint and store them in the flattened array
        jacobian_data[jacobian_index] = bar1s[i].jacobian_entry_x(mass_list);
        jacobian_data[jacobian_index + 1] = bar1s[i].jacobian_entry_y(mass_list);
    }

    // Loop through each rigid_bar_2 constraint
    for (int i = 0; i < num_bar2s; i++)
    {
        // Calculate the index in the Matrix_Block list for the current constraint
        int bar2Ind = num_bar1s + i;

        // Calculate the index in the flattened Jacobian array for the current constraint
        int jacobian_index = 2 * num_bar1s + 4 * i;

        // Create a Matrix_Block for the current constraint and store its information
        jacobian[bar2Ind] = Matrix_Block(1, 4, bar2Ind, 2 * bar2s[i].attached_masses[0], jacobian_index);

        // Compute the Jacobian entries for the current constraint and store them in the flattened array
        jacobian_data[jacobian_index] = bar2s[i].jacobian_entry_x1(mass_list);
        jacobian_data[jacobian_index + 1] = bar2s[i].jacobian_entry_y1(mass_list);
        jacobian_data[jacobian_index + 2] = bar2s[i].jacobian_entry_x2(mass_list);
        jacobian_data[jacobian_index + 3] = bar2s[i].jacobian_entry_y2(mass_list);
    }
}


// WORKING CREATE_JACOBIAN
//void State_Getter::create_jacobian_old(Circular_Rigid_Body *mass_list, Rigid_Bar_1 *bar1s, int num_bar1s, Rigid_Bar_2 *bar2s, int num_bar2s)
//{
//	// Jacobian is made up of a list of matrix blocks representing nonzero values. There will be as many matrix blocks as there are constraints
//	// Dynamically allocate space for the matrix block list m x 2n
//	jacobian = (Matrix_Block *)malloc(num_constraints * sizeof(Matrix_Block));
//
//	// Now I need to create my matrix blocks
//	// Loop through each list of constraints
//
//	// Loop through lists of each type of constraint, these are all the constraint functions
//	// Rigid_Bar_1s.
//	for (int i = 0; i < num_bar1s; i++)
//	{
//		// Each rigid_bar_1 only applies a partial derivative to the x and y of 1 attached mass. this contributes two nonzero terms to the jacobian
//		jacobian[i].rows = 1;
//		jacobian[i].cols = 2;
//		// Row offset: 1 row per constraint object
//		jacobian[i].row = i;
//		// Column offset: column of origin of matrix block starts at the mass attached to the rigid rod.
//		// 2 derivatives per mass so multiplied by 2.
//		jacobian[i].col = 2 * bar1s[i].attached_mass;
//
//		// Allocate the space needed for the matrix in the matrix block.
//		// Rigid_Bar_1 constraints contribute 2 non zero terms for the 1 particle affected. partial x and partial y
//		jacobian[i].matrix = (double *)malloc(2 * sizeof(double));
//
//		// Fill up the matrix
//		//jacobian[count].matrix[0] = bar1s[i].partial_x(mass_list);
//		jacobian[i].matrix[0] = bar1s[i].jacobian_entry_x(mass_list);
//		jacobian[i].matrix[1] = bar1s[i].jacobian_entry_y(mass_list);
//	}
//	// Rigid_Bar_2s
//	if (num_bar2s > 0)
//	{
//		for (int i = 0; i < num_bar2s; i++)
//		{
//			int bar2Ind = num_bar1s + i;
//			// Each rigid_bar_2 applies a partial derivative to the x and y of both attached masses. this contributes four nonzero terms to the jacobian
//			jacobian[bar2Ind].rows = 1;
//			jacobian[bar2Ind].cols = 4;
//			// Row offset: 1 row per constraint object
//			jacobian[bar2Ind].row = num_bar1s + i;
//			// Column offset: column of origin of matrix block starts at the first mass the rigid rod is attached to.
//			// 2 derivatives per mass so multiplied by 2.
//			// Assume minimum attached mass is first in the attached mass list
//			jacobian[bar2Ind].col = 2 * bar2s[i].attached_masses[0];
//			// Allocate the space needed for the matrix in the matrix block.
//			// Rigid_Bar_2 constraints contribute 4 non zero terms for the 2 particles affected. Need 4 total terms. x1, y1, x2, y2.
//			jacobian[bar2Ind].matrix = (double *)malloc(4 * sizeof(double));
//
//			// Fill up the matrix
//			jacobian[bar2Ind].matrix[0] = bar2s[i].jacobian_entry_x1(mass_list);
//			jacobian[bar2Ind].matrix[1] = bar2s[i].jacobian_entry_y1(mass_list);
//			jacobian[bar2Ind].matrix[2] = bar2s[i].jacobian_entry_x2(mass_list);
//			jacobian[bar2Ind].matrix[3] = bar2s[i].jacobian_entry_y2(mass_list);
//		}
//	}
//}

// Testing new create jacobian derivative with new matrix block style: NEEDS TESTING
void State_Getter::create_jacobian_derivative(Circular_Rigid_Body *mass_list, Rigid_Bar_1 *bar1s, int num_bar1s, Rigid_Bar_2 *bar2s, int num_bar2s)
{
    // Allocate memory for the matrix block list, one block for each constraint
    jacobian_derivative = (Matrix_Block *)malloc(num_constraints * sizeof(Matrix_Block));

    // Allocate memory for the entire Jacobian derivative array at once
    jacobian_derivative_data = (double *)malloc((2 * num_bar1s + num_bar2s * 4) * sizeof(double));

    // Loop through Rigid_Bar_1 constraints
    for (int i = 0; i < num_bar1s; i++)
    {
        // Calculate the index in the flattened Jacobian derivative array for the current constraint
        int jacobian_derivative_index = 2 * i;

        // Create a Matrix_Block for the current constraint and store its information
        jacobian_derivative[i] = Matrix_Block(1, 2, i, 2 * bar1s[i].attached_mass, jacobian_derivative_index);

        // Compute the Jacobian derivative entries for the current constraint and store them in the flattened array
        jacobian_derivative_data[jacobian_derivative_index] = bar1s[i].jacobian_derivative_entry_x(mass_list);
        jacobian_derivative_data[jacobian_derivative_index + 1] = bar1s[i].jacobian_derivative_entry_y(mass_list);
    }

    // Loop through Rigid_Bar_2 constraints
    for (int i = 0; i < num_bar2s; i++)
    {
        // Calculate the index in the Matrix_Block list for the current constraint
        int bar2Ind = num_bar1s + i;

        // Calculate the index in the flattened Jacobian derivative array for the current constraint
        int jacobian_derivative_index = 2 * num_bar1s + 4 * i;

        // Create a Matrix_Block for the current constraint and store its information
        jacobian_derivative[bar2Ind] = Matrix_Block(1, 4, bar2Ind, 2 * bar2s[i].attached_masses[0], jacobian_derivative_index);

        // Compute the Jacobian derivative entries for the current constraint and store them in the flattened array
        jacobian_derivative_data[jacobian_derivative_index] = bar2s[i].jacobian_derivative_entry_x1(mass_list);
        jacobian_derivative_data[jacobian_derivative_index + 1] = bar2s[i].jacobian_derivative_entry_y1(mass_list);
        jacobian_derivative_data[jacobian_derivative_index + 2] = bar2s[i].jacobian_derivative_entry_x2(mass_list);
        jacobian_derivative_data[jacobian_derivative_index + 3] = bar2s[i].jacobian_derivative_entry_y2(mass_list);
    }
}

// Working create_jacobian_derivative
//void State_Getter::create_jacobian_derivative_old(Circular_Rigid_Body *mass_list, Rigid_Bar_1 *bar1s, int num_bar1s, Rigid_Bar_2 *bar2s, int num_bar2s)
//{
//	// Jacobian is made up of a list of matrix blocks representing nonzero values. There will be as many matrix blocks as there are constraints
//	// Dynamically allocate space for the matrix block list m x 2n
//	jacobian_derivative = (Matrix_Block *)malloc(num_constraints * sizeof(Matrix_Block));
//
//	// Now I need to create my matrix blocks
//	// Loop through each list of constraints
//	// Loop through lists of each type of constraint, these are all the constraint functions
//	// Rigid_Bar_1s.
//	for (int i = 0; i < num_bar1s; i++)
//	{
//		// Each rigid_bar_1 only applies a double partial derivative multiplied by velocity term for both x and y of 1 attached mass.
//		// this contributes two nonzero terms to the jacobian
//		jacobian_derivative[i].rows = 1;
//		jacobian_derivative[i].cols = 2;
//		// Row offset: 1 row per constraint object
//		jacobian_derivative[i].row = i;
//		// Column offset: column of origin of matrix block starts at the mass attached to the rigid rod.
//		// 2 derivatives per mass so multiplied by 2.
//		jacobian_derivative[i].col = 2 * bar1s[i].attached_mass;
//
//		// Allocate the space needed for the matrix in the matrix block.
//		// Rigid_Bar_1 constraints contribute 2 non zero terms for the 1 particle affected. partial x and partial y
//		jacobian_derivative[i].matrix = (double *)malloc(2 * sizeof(double));
//
//		// Fill up the matrix
//		//jacobian[count].matrix[0] = bar1s[i].partial_x(mass_list);
//		jacobian_derivative[i].matrix[0] = bar1s[i].jacobian_derivative_entry_x(mass_list);
//		jacobian_derivative[i].matrix[1] = bar1s[i].jacobian_derivative_entry_y(mass_list);
//	}
//	if (num_bar2s > 0)
//	{
//		// Rigid_Bar_2s
//		for (int i = 0; i < num_bar2s; i++)
//		{
//			// Each rigid_bar_2 applies a partial derivative to the x and y of both attached masses. this contributes four nonzero terms to the jacobian
//			jacobian_derivative[num_bar1s + i].rows = 1;
//			jacobian_derivative[num_bar1s + i].cols = 4;
//			// Row offset: 1 row per constraint object
//			jacobian_derivative[num_bar1s + i].row = num_bar1s + i;
//			// Column offset: column of origin of matrix block starts at the first mass the rigid rod is attached to.
//			// 2 derivatives per mass so multiplied by 2.
//			// Assume minimum attached mass is first in the attached mass list
//			jacobian_derivative[num_bar1s + i].col = 2 * bar2s[i].attached_masses[0];
//			// Allocate the space needed for the matrix in the matrix block.
//			// Rigid_Bar_2 constraints contribute 4 non zero terms for the 2 particles affected. Need 4 total terms. x1, y1, theta1, x2, y2.
//			jacobian_derivative[num_bar1s + i].matrix = (double *)malloc(4 * sizeof(double));
//
//			// Fill up the matrix
//			jacobian_derivative[num_bar1s + i].matrix[0] = bar2s[i].jacobian_derivative_entry_x1(mass_list);
//			jacobian_derivative[num_bar1s + i].matrix[1] = bar2s[i].jacobian_derivative_entry_y1(mass_list);
//			jacobian_derivative[num_bar1s + i].matrix[2] = bar2s[i].jacobian_derivative_entry_x2(mass_list);
//			jacobian_derivative[num_bar1s + i].matrix[3] = bar2s[i].jacobian_derivative_entry_y2(mass_list);
//		}
//	}
//}


// New implementation for new Matrix Block class
// Calculate the B vector evaluating the right side: -J_dot * q_dot - J * M^-1 * F_ext
void State_Getter::calculate_b(double *b)
{
	// First calculate the left side. J_dot * q_dot
	double *lhs = new double[num_constraints]();

	matrix_blocks_times_vector(jacobian_derivative, num_constraints, jacobian_derivative_data, state_vector_derivative, lhs);

	// Then multiply it by negative 1.
	for (int i = 0; i < num_constraints; i++)
	{
		lhs[i] = -1 * lhs[i];
	}

	// Now do the right hand side.
	// First inverse mass times F_ext
	double *first_res = new double[2 * num_bodies]();
	diagonal_times_vector(inverse_mass_matrix, force_ext_vector, 2 * num_bodies, first_res);

	// Then J times that
	double *rhs = new double[num_constraints]();
	matrix_blocks_times_vector(jacobian, num_constraints, jacobian_data, first_res, rhs);

	// b = lhs - rhs
	for (int i = 0; i < num_constraints; i++)
	{
		b[i] = lhs[i] - rhs[i];
	}

	// Free stuff
	delete[] lhs;
	delete[] rhs;
	delete[] first_res;
}

// Working implementation with old Matrix Block class
// Calculate the B vector evaluating the right side: -J_dot * q_dot - J * M^-1 * F_ext
//void State_Getter::calculate_b(double *b)
//{
//	// First calculate the left side. J_dot * q_dot
//	double *lhs = new double[num_constraints]();
//
//	matrix_blocks_times_vector(jacobian_derivative, num_constraints, state_vector_derivative, lhs);
//
//	// Then multiply it by negative 1.
//	for (int i = 0; i < num_constraints; i++)
//	{
//		lhs[i] = -1 * lhs[i];
//	}
//
//	// Now do the right hand side.
//	// First inverse mass times F_ext
//	double *first_res = new double[2 * num_bodies]();
//	diagonal_times_vector(inverse_mass_matrix, force_ext_vector, 2 * num_bodies, first_res);
//
//	// Then J times that
//	double *rhs = new double[num_constraints]();
//	matrix_blocks_times_vector(jacobian, num_constraints, first_res, rhs);
//
//	// b = lhs - rhs
//	for (int i = 0; i < num_constraints; i++)
//	{
//		b[i] = lhs[i] - rhs[i];
//	}
//
//	// Free stuff
//	delete[] lhs;
//	delete[] rhs;
//	delete[] first_res;
//}

void State_Getter::calculate_net_force_vector(void)
{
	// Allocate the vector pointer
	net_force_vector = (double *)malloc(2 * num_bodies * sizeof(double));

	for (int i = 0; i < 2 * num_bodies; i++)
	{
		net_force_vector[i] = force_ext_vector[i] + constraint_force_vector[i];
	}
}

