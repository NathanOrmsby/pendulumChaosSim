/*
 * toFile.cpp
 *
 *  Created on: Feb 21, 2023
 *      Author: norms
 */

#include "../headers/toFile.h"
#include "../headers/utils.h"
#include "../headers/DoublePendulum.h"
#include <vector>
#include <cmath>
#include <tuple>
#include <iostream>
#include <fstream>
#include <string>


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

// Write the grid of angles to csv in degrees
void writeAnglesToCSV(double *angles, int resolution[2])
{
    // Open the output file
	std::string fpath = "D:\\pendulumSimulator\\pendulumData\\";
	std::string filename = "angleGrid.csv";
    std::ofstream outFile(fpath + filename);

    // Write the header row
    outFile << "theta1,theta2" << std::endl;

    // Write the angle data
    for (int i = 0; i < resolution[0] * resolution[1]; i++)
    {
        // Write the theta1 angle
        outFile << radiansToDegrees(angles[2 * i]) << ",";

        // Write the theta2 angle
        outFile << radiansToDegrees(angles[2 * i + 1]) << std::endl;
    }

    // Close the output file
    outFile.close();
}

// Plot the initial angles of all pendulums and their perturbations in the simulation

std::vector<std::vector<double>> DoublePendulumAnglesFromMassListAndWriteToCSV(Circular_Rigid_Body *mass_list, int mass_list_size) {

	std::string fpath = "D:\\pendulumSimulator\\pendulumData\\";
	std::string filename = "lyapunovPendulumSystems.csv";
	std::ofstream outfile(fpath + filename);
	std::vector<std::vector<double>> pendulum_angles;

    if (mass_list_size % 2 != 0) {
        // Invalid mass list size
        return pendulum_angles;
    }


    if (!outfile.is_open()) {
        std::cerr << "Error: Unable to open file for writing: " << filename << std::endl;
        return pendulum_angles;
    }

    for (int i = 0; i < mass_list_size; i += 2) {
        Circular_Rigid_Body &mass1 = mass_list[i];
        Circular_Rigid_Body &mass2 = mass_list[i + 1];
        std::vector<double> angles = DoublePendulumAnglesFromMasses(mass1, mass2);
        angles.front() = radiansToDegrees(angles.front());
        angles.back() = radiansToDegrees(angles.back());
        pendulum_angles.push_back(angles);

        if (angles.size() >= 2) {
            outfile << angles[0] << "," << angles[1] << "\n";
        }
    }

    outfile.close();
    return pendulum_angles;
}



