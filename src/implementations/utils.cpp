#include <iostream>

#include "../headers/utils.h"
#include <cmath>
int minimum(int a, int b)
{
	if (a < b)
	{
		return a;
	}
	else if (a > b)
	{
		return b;
	}
	// If they're equal
	else
	{
		return 9999;
	}
}


int maximum(int a, int b)
{
	if (a > b)
	{
		return a;
	}
	else if (a < b)
	{
		return b;
	}
	// If they're equal
	else
	{
		return 9999;
	}
}

double degreesToRadians(double degrees)
{
    return degrees * M_PI / 180.0;
}

double radiansToDegrees(double radians) {
	return radians * (180 / M_PI);
}

// Not in use with new style of Matrix Block
// ---------------------------------------------------------------
//void free_matrix_block_list(Matrix_Block **block_list, int num_blocks)
//{
//	for (int i = 0; i < num_blocks; i++)
//	{
//		// Free the matrix arrays in each block
//		free(block_list[i]->matrix);
//	}
//
//	free(block_list);
//}
//
//void print_matrix_block_list(Matrix_Block *block_list, int num_blocks)
//{
//	for (int i = 0; i < num_blocks; i++)
//	{
//		std::cout << "Matrix Block " << i << std::endl;
//		std::cout << "Offset: row: " << block_list[i].row << " col: " << block_list[i].col << std::endl;
//		std::cout << "Printing matrix: " << std::endl;
//		for (int j = 0; j < block_list[i].cols; j++)
//		{
//			std::cout << block_list[i].matrix[j] << " ";
//		}
//		std::cout << std::endl;
//	}
//}
