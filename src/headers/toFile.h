/*
 * toFile.h
 *
 *  Created on: Feb 21, 2023
 *      Author: norms
 */

#ifndef HEADERS_TOFILE_H_
#define HEADERS_TOFILE_H_

#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <sstream>

#include "vectors.h"

// Writes the pendulum position data to file. The pivot position, and positions of all masses
void pendulumToFile(Vector pivotData, Vector **massData, int dataLen, int num_masses);


#endif /* HEADERS_TOFILE_H_ */
