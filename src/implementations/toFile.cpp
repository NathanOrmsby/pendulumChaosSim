/*
 * toFile.cpp
 *
 *  Created on: Feb 21, 2023
 *      Author: norms
 */

#include "../headers/toFile.h"

// Writes the pendulum position data to file. The pivot position, and positions of all masses
void pendulumToFile(Vector pivotData, Vector **massData, int dataLen, int num_masses)
{
	// File path
	std::string fpath = "D:\\pendulumSimulator\\pendulumData\\";

	// File name
	std::string f = "pendulum.csv";

	// Initialize the output file streams
	std::ofstream file;

	// Open the file
	file.open(fpath + f);

	// File header

	// Pivot
	file << "pivotx,pivoty,";

	// Masses
	for (int i = 0; i < num_masses - 1; ++i)
	{
		file << "mx" + std::to_string(i) + "," + "my" + std::to_string(i) + ",";
	}
	file << "mx" + std::to_string(num_masses - 1) + "," + "my" + std::to_string(num_masses - 1) + "\n";

	// Add data to the file
	for (int i = 0; i < dataLen; ++i)
	{
		// Pivot
		file << std::to_string(pivotData.x) << "," << std::to_string(pivotData.y) << ",";

		// Masses
		for (int j = 0; j < num_masses - 1; ++j)
		{
			file << std::to_string(massData[i][j].x) + "," + std::to_string(massData[i][j].y) + ",";
		}
		file << std::to_string(massData[i][num_masses - 1].x) + "," + std::to_string(massData[i][num_masses - 1].y) + "\n";
	}

	// Close the file
	file.close();
}


