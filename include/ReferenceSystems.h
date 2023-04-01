#ifndef REFERENCESYSTEMS_H
#define REFERENCESYSTEMS_H

#include "LinearAlgebra.h"

// ----------------------------------------------
// GCRF2Perifocal
// ----------------------------------------------
// Purpose:
//    Transformation matrix from GCRF (Geocentric Celestial Reference Frame) to Perifocal reference frame (z axis peripendicular to the orbit plane and x axis in the perigee direction)
//
// Parameters:
//    i: inclination of the orbit [rad]
//    Omega: right ascension of the ascending node [rad]
//    omega: argument of perigee [rad]
//
// Returns:
//    Transformation matrix from ITRF to Perifocal
// ----------------------------------------------
Matrix GCRF2Perifocal(double i, double Omega, double omega);

// ----------------------------------------------
// GCRF2ITRF
// ----------------------------------------------
// Purpose:
//    Transformation matrix from GCRF (Geocentric Celestial Reference Frame) to ITRF (International Terrestrial Reference Frame)
//
// Parameters:
//    mjd_TT: time in Modified Julian date TT (terrestrial time)
//
// Returns:
//    Transformation matrix from ICRF to ITRF
// ----------------------------------------------
Matrix GCRF2ITRF(double mjd_TT);

// ----------------------------------------------
// MeanObliquity
// ----------------------------------------------
// Purpose:
//    Calculate mean obliquity of the ecliptic
//
// Parameters:
//    mjd_TT: time in Modified Julian date TT (terrestrial time)
//
// Returns:
//    Mean obliquity of the ecliptic [rad]
// ----------------------------------------------
double MeanObliquity(double mjd_TT);

// ----------------------------------------------
// EqEquinoxes
// ----------------------------------------------
// Purpose:
//    Calculate equation of the equinoxes
//
// Parameters:
//    mjd_TT: time in Modified Julian date TT (terrestrial time)
//
// Returns:
//    Equation of the equinoxes [rad]
// ----------------------------------------------
double EqEquinoxes(double mjd_TT);

// ----------------------------------------------
// GMST
// ----------------------------------------------
// Purpose:
//    Calculate Greenwich mean sidereal time
//
// Parameters:
//    mjd_UT1: time in Modified Julian date UT1 (universal time)
//
// Returns:
//    Greenwich mean sidereal time [rad]
// ----------------------------------------------
double GMST(double mjd_UT1);

// ----------------------------------------------
// GAST
// ----------------------------------------------
// Purpose:
//    Calculate Greenwich apparent sidereal time
//
// Parameters:
//    mjd_UT1: time in Modified Julian date UT1 (universal time)
//
// Returns:
//    Greenwich apparent sidereal time [rad]
// ----------------------------------------------
double GAST(double mjd_UT1);

// ----------------------------------------------
// get_gregorian_date
// ----------------------------------------------
// Purpose:
//    Convert TLE format date to Gregorian date
//
// Parameters:
//    year: year in YYYY format
//    day_frac: day fraction of the year in TLE format
//    &month: where to store month
//    &day: where to store day
//
// Returns:
//    None
// ----------------------------------------------
void get_gregorian_date(int year, double day_frac, int& month, double& day);

// ----------------------------------------------
// get_mjd_TT
// ----------------------------------------------
// Purpose:
//    Convert Gregorian date to Modified Julian date TT (terrestrial time)
//
// Parameters:
//    yyyy: year in YYYY format
//    mm: month
//    dd: day
//
// Returns:
//    Modified Julian date TT (terrestrial time)
// ----------------------------------------------
double get_mjd_TT(double yyyy, double mm, double dd);

// ----------------------------------------------
// PrecessionMatrix
// ----------------------------------------------
// Purpose:
//    Calculate precession matrix that aligns the mean equator and mean vernal equinox of the J2000.0 reference frame with the mean equator and mean vernal equinox of the given time
//
// Parameters:
//    mjd_TT: time in Modified Julian date TT (terrestrial time)
//
// Returns:
//    Precession matrix (a product of three rotation matrices)
// ----------------------------------------------
Matrix PrecessionMatrix(double mjd_TT);

// ----------------------------------------------
// NutationAngles
// ----------------------------------------------
// Purpose:
//    Calculate nutation angles dpsi (nutation in longitude) and deps (nutation in obliquity)
//
// Parameters:
//    mjd_TT: time in Modified Julian date TT (terrestrial time)
//    &dpsi: where to store nutation in longitude [rad]
//    &deps: where to store nutation in obliquity [rad]
//
// Returns:
//    None
// ----------------------------------------------

void NutationAngles(double mjd_TT, double& dpsi, double& deps);

// ----------------------------------------------
// NutationMatrix
// ----------------------------------------------
// Purpose:
//    Calculate nutation matrix that aligns the mean equator and mean vernal equinox of the given time with the true equator and true vernal equinox of the given time
//
// Parameters:
//    mjd_TT: time in Modified Julian date TT (terrestrial time)
//
// Returns:
//    Nutation matrix (a product of three rotation matrices)
// ----------------------------------------------
Matrix NutationMatrix(double mjd_TT);

// ----------------------------------------------
// RotationMatrix
// ----------------------------------------------
// Purpose:
//    Calculate rotation matrix that aligns the true vernal equinox with th 'mean' (up to polar motion) prime meridian of the given time
//
// Parameters:
//    mjd_UT1: time in Modified Julian date UT1 (universal time)
//
// Returns:
//    Rotation matrix
// ----------------------------------------------
Matrix RotationMatrix(double mjd_UT1);

// ----------------------------------------------
// PolarMotion
// ----------------------------------------------
// Purpose:
//    Calculate polar motion angles xp and yp
//
// Parameters:
//    mjd_TT: time in Modified Julian date TT (terrestrial time)
//    &xp: where to store polar motion angle xp [rad]
//    &yp: where to store polar motion angle yp [rad]
//
// Returns:
//    None
// ----------------------------------------------
void PolarMotion(double mjd_TT, double& xp, double& yp);

// ----------------------------------------------
// PolarMotionMatrix
// ----------------------------------------------
// Purpose:
//    Calculate polar motion matrix that corrects the errors in the Earth's rotation axis
//
// Parameters:
//    mjd_TT: time in Modified Julian date TT (terrestrial time)
//
// Returns:
//    Polar motion matrix (a product of two rotation matrices)
// ----------------------------------------------
Matrix PolarMotionMatrix(double mjd_TT);

#endif  // REFERENCESYSTEMS_H