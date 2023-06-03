#ifndef REFERENCESYSTEMS_H
#define REFERENCESYSTEMS_H

#include "LinearAlgebra.h"

// ----------------------------------------------
// Perifocal2ECI
// ----------------------------------------------
// Purpose:
//    Transformation matrix from Perifocal reference frame (z axis peripendicular to the orbit plane and x axis in the perigee direction) to ECI (Earth-centered interial frame) of J2000
//
// Parameters:
//    i: inclination of the orbit [rad]
//    Omega: right ascension of the ascending node [rad]
//    omega: argument of perigee [rad]
//
// Returns:
//    Transformation matrix from Perifocal to ECI
// ----------------------------------------------
Matrix Perifocal2ECI(double i, double Omega, double omega);

// ----------------------------------------------
// J20002ECEF
// ----------------------------------------------
// Purpose:
//    Transformation matrix from ECI (Earth-centered inertial) of J2000 to ECEF (Earth-centered Earth-fixed)
//
// Parameters:
//    mjd_TT: time in Modified Julian date TT (terrestrial time)
//
// Returns:
//    Transformation matrix from ECI to ECEF
// ----------------------------------------------
Matrix J20002ECEF(double mjd_TT);

// ----------------------------------------------
// ECI2TEME
// ----------------------------------------------
// Purpose:
//    Transformation matrix from TEME (True Equator Mean Equinox) to ECI (Earth-centered inertial) of J2000
//
// Parameters:
//    mjd_TT: time in Modified Julian date TT (terrestrial time)
//
// Returns:
//    Transformation matrix from TEME to ECI
// ----------------------------------------------
Matrix ECI2TEME(double mjd_TT);

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
//    mjd_TT: time in Modified Julian date TT (terrestrial time)
//
// Returns:
//    Greenwich apparent sidereal time [rad]
// ----------------------------------------------
double GAST(double mjd_TT);

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
// get_mjd
// ----------------------------------------------
// Purpose:
//    Convert Gregorian date (in UT1 or TT) to Modified Julian date (in UT1 or TT, respectively)
//
// Parameters:
//    yyyy: year in yyyy format
//    mm: month in mm format
//    dd: day in dd format
//
// Returns:
//    Modified Julian date TT (terrestrial time)
// ----------------------------------------------
double get_mjd(double yyyy, double mm, double dd);

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
// EquinoxMatrix
// ----------------------------------------------
// Purpose:
//    Calculate equinox matrix that aligns the mean equator with the true equator of the given time
//
// Parameters:
//    mjd_TT: time in Modified Julian date TT (terrestrial time)
//
// Returns:
//    Equinox matrix
// ----------------------------------------------
Matrix EquinoxMatrix(double mjd_TT);

// ----------------------------------------------
// RotationMatrix
// ----------------------------------------------
// Purpose:
//    Calculate rotation matrix that aligns the true vernal equinox with the 'mean' (up to polar motion) prime meridian of the given time
//
// Parameters:
//    mjd_TT: time in Modified Julian date TT (terrestrial time)
//
// Returns:
//    Rotation matrix
// ----------------------------------------------
Matrix RotationMatrix(double mjd_TT);

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

// ----------------------------------------------
// Sun
// ----------------------------------------------
// Purpose:
//    Calculate position vector of the Sun in the J2000.0 reference frame
//
// Parameters:
//    Mjd_TT: time in Modified Julian date TT (terrestrial time)
//
// Returns:
//    Position vector of the Sun in the J2000.0 reference frame [m]
// ----------------------------------------------
Vector Sun(double Mjd_TT);

// ----------------------------------------------
// Moon
// ----------------------------------------------
// Purpose:
//    Calculate position vector of the Moon in the J2000.0 reference frame
//
// Parameters:
//    Mjd_TT: time in Modified Julian date TT (terrestrial time)
//
// Returns:
//    Position vector of the Moon in the J2000.0 reference frame [m]
// ----------------------------------------------
Vector Moon(double Mjd_TT);

#endif  // REFERENCESYSTEMS_H