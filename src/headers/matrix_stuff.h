/*
 * matrix_stuff.h
 *
 *  Created on: Nov 8, 2022
 *      Author: norms
 */




#ifndef MATRIX_STUFF_H_
#define MATRIX_STUFF_H_

#include "get_state.h"

class State_Getter;

// Optimized Matrix Block: Operates on an external matrix double array
class Matrix_Block {
public:
    int rows;
    int cols;
    int row;
    int col;
    int jacobian_index; // index into the jacobian array

    // Constructor
    Matrix_Block(int r, int c, int row_offset, int col_offset, int jacobian_index) : rows(r), cols(c), row(row_offset), col(col_offset) , jacobian_index(jacobian_index) {}
};

// Matrix block, represents only non zero values and saves time when matrices have a lot of zeros
//class Matrix_Block
//{
//	public:
//	// row and column number of origin. this is the offset
//	int row;
//	int col;
//	// row and column lengths.
//	// jlength is always the number of dimensions the particle is in. in this case, 3. ilength is however many constraints that are of the same type we have. Number of rows.
//	int rows;
//	int cols;
//	// double array containing the matrix. matrix will be ilength * jlength. Number of columns of matrix. 1d array.
//	double *matrix;
//};

// -------------------------------------------- MATRIX FUNCTIONS USED BY ... --------------------------------------------


//-------------------------------------------- Used by BICONJUGATE GRADIENT --------------------------------------------

// Function that calculates A*x: A is always symmetric in our problem. A= JM^-1JT
void A_times_x(State_Getter *state, double *x, double *result);

// Function that either returns the magnitude of the vector, or the magnitude of the largest component.
// itol is the protocol
double calculate_relevant_norm(double *v, int n, int itol);

// Function that solves the Ax = B for current minimizers
void solve_for_x(int n, double *A_diagonal, double *b, double *x);

//-------------------------------------------- Used by STATE GETTER --------------------------------------------

// Multiplies a list of matrix blocks times a vector. Results in a vector
void matrix_blocks_times_vector(Matrix_Block *block_list, int num_blocks, double *jacobian, double *vector, double *result);

// Diagonal matrix multiplied by vector
void diagonal_times_vector(double *diagonal, double *vector, int n, double *result);

// Multiplies the transpose of a list of matrix blocks by a vector. Results in a 1d vector
void matrix_blocks_transpose_times_vector(Matrix_Block *block_list, int num_blocks, double *jacobian, double *vector, double *result);

// Transpose matrix blocks
//void matrix_blocks_transpose(Matrix_Block *block_list, int num_blocks, Matrix_Block *result_blocks);

// NOT USED

//// Returns the transpose of a symmetric matrix
//void matrix_transpose_symmetric(double *A, int rows, int cols, double *transposed_matrix);
//
//// Returns the transpose of a nonsymmetric matrix
//void matrix_transpose_nonsymmetric(double *A, int rows, int cols, double *matrix_transposed);
//
//// Matrix multiplication
//void matrix_multiplication_ROW_BY_ROW(double *A, double *B, int *A_dim, int *B_dim, double *result);
//
//// Normal matrix multiplication
//void matrix_multiplication(double *A, double *B, int *A_dim, int *B_dim, double *result);
//
//// Take inverse of a diagonal matrix, represented as a vector
//void inverse_diagonal_matrix(double *matrix, int length, double *result);
//
//// Jacobian times the diagonal inverse mass matrix. Edits a list of resulting matrix blocks as the resulting m x 3n matrix is sparse.
//void matrix_blocks_times_diagonal(Matrix_Block *block_list, int num_blocks, double *inverse_mass_matrix, Matrix_Block *result_blocks);
//
//
//

//
//
//
//

//
// zeros all elements in array
void zero_vector(double  *vector, int length);

// Fills array with ones
void ones_vector(double *vector, int length);

#endif /* MATRIX_STUFF_H_ */
