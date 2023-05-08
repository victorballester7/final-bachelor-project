#ifndef MISC_H
#define MISC_H

#include <fstream>

#include "constants.h"

// ----------------------------------------------
// to2Pi
// ----------------------------------------------
// Purpose:
//    Transform angle to the interval [0, 2*pi)
//
// Parameters:
//    angle: angle to be transformed
//
// Returns:
//    Transformed angle
// ----------------------------------------------
double to2Pi(double angle);
int numLines(std::string filename);
#endif  // MISC_H