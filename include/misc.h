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

// ----------------------------------------------
// numLines
// ----------------------------------------------
// Purpose:
//    Count number of lines in a file
//
// Parameters:
//    filename: name of the file
//
// Returns:
//    Number of lines in the file
// ----------------------------------------------
int numLines(std::string filename);

// ----------------------------------------------
// mean
// ----------------------------------------------
// Purpose:
//    Compute the mean of an array
//
// Parameters:
//    array: array of doubles
//    n: number of elements in the array
//
// Returns:
//    Mean of the array
// ----------------------------------------------
double mean(double *array, int n);

// ----------------------------------------------
// variance
// ----------------------------------------------
// Purpose:
//    Compute the variance of an array
//
// Parameters:
//    array: array of doubles
//    n: number of elements in the array
//
// Returns:
//    Variance of the array
// ----------------------------------------------
double variance(double *array, int n);
#endif  // MISC_H