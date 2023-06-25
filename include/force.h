#ifndef FORCE_H
#define FORCE_H

#include "LinearAlgebra.h"
#include "Satellite.h"

typedef struct {
  int n_max;          // maximum degree of the spherical harmonics expansion(0...N_JDM3 - 1)
  int m_max;          // maximum order of the spherical harmonics expansion (0...n_max). m_max = 0 for the zonal harmonics only
  bool pointEarth;    // true if point mass approximation of the Earth is used, false if spherical harmonics expansion is used
  bool sun;           // true if the Sun is included, false otherwise
  bool moon;          // true if the Moon is included, false otherwise
  bool otherPlanets;  // true if the other planets are included, false otherwise
  bool solar_rad;     // true if solar radiation pressure is included, false otherwise
  bool atmo_drag;     // true if solar radiation pressure is included, false otherwise
  double Am;          // Area-to-mass ratio of the satellite: A/m [m^2 / kg] cross-sectional area of the satellite [m^2]
} args_gravField;

// ----------------------------------------------
// AccelPointEarth
// ----------------------------------------------
// Purpose:
//    Calculate acceleration of the satellite at position r due the gravity exerced by the Earth considered as a point mass
//
// Parameters:
//    r: position vector of the satellite in the inertial reference frame [m]
//
// Return value:
//    acceleration vector in the inertial reference frame [m/s^2]
// ----------------------------------------------
Vector AccelPointEarth(const Vector& r);

// ----------------------------------------------
// AccelPointMass
// ----------------------------------------------
// Purpose:
//    Calculate acceleration of the satellite at position r due the gravity exerced by a point mass located at s
//
// Parameters:
//    r: position vector of the satellite in the inertial reference frame [m]
//    s: position vector of the point mass in the inertial reference frame [m]
//    GM: gravitational parameter of the central body [m^3/s^2]
//
// Return value:
//    acceleration vector in the inertial reference frame [m/s^2]
// ----------------------------------------------
Vector AccelPointMass(const Vector& r, const Vector& s, double GM);

// ----------------------------------------------
// AccelSolarRad
// ----------------------------------------------
// Purpose:
//    Calculate acceleration of the satellite at position r due the solar radiation pressure
//
// Parameters:
//    r: position vector of the satellite in the inertial reference frame [m]
//    r_Sun: position vector of the Sun in the inertial reference frame [m]
//    Area_mass: Area-to-mass ratio of the satellite: A/m [m^2 / kg]
//
// Return value:
//    acceleration vector in the inertial reference frame [m/s^2]
// ----------------------------------------------
Vector AccelSolarRad(const Vector& r, const Vector& r_Sun, double Area_mass);

// ----------------------------------------------
// Illumination
// ----------------------------------------------
// Purpose:
//    Calculate the fraction of the satellite's surface illuminated by the Sun
//
// Parameters:
//    r: position vector of the satellite in the inertial reference frame [m]
//    r_Sun: position vector of the Sun in the inertial reference frame [m]
//
// Return value:
//   0: satellite is in the shadow of the Earth
//   1: satellite is fully illuminated by the Sun
double Illumination(const Vector& r, const Vector& r_Sun);

// ----------------------------------------------
// AccelDrag
// ----------------------------------------------
// Purpose:
//    Calculate acceleration of the satellite at position r due the atmospheric drag
//
// Parameters:
//    r: position vector of the satellite in the inertial reference frame of J2000 [m]
//    v: velocity vector of the satellite in the inertial reference frame of J2000 [m/s]
//    mjd_TT: modified Julian date (Terrestrial Time)
//    NP: transformation matrix from the inertial reference frame of J2000 to the reference frame true of date (without the precession and nutation matrices)
//    Area_mass: Area-to-mass ratio of the satellite: A/m [m^2 / kg]
//
// Return value:
//    acceleration vector in the inertial reference frame of J2000 [m/s^2]
// ----------------------------------------------
Vector AccelDrag(const Vector& r, const Vector& v, double mjd_TT, const Matrix& NP, double Area_mass);

// ----------------------------------------------
// Density_HP
// ----------------------------------------------
// Purpose:
//    Calculate the atmospheric density at the altitude of the satellite using the Harris-Priester model
//
// Parameters:
//    mjd_TT: modified Julian date (Terrestrial Time)
//    r_NP: position vector of the in the reference frame true of date (without the precession and nutation matrices) [m]
//
// Return value:
//    atmospheric density [kg/m^3]
// ----------------------------------------------
double Density_HP(double mjd_TT, const Vector& r_NP);

// ----------------------------------------------
// AccelHarmonic
// ----------------------------------------------
// Purpose:
//    Calculate acceleration due to gravity using the spherical harmonics expansion of the gravitational potential
//
// Parameters:
//    r: position vector in the inertial reference frame [m]
//    E: transformation matrix from the inertial reference frame to the reference frame of the spherical harmonics expansion (fixed to the Earth)
//    n_max: maximum degree of the spherical harmonics expansion (0...N_JDM3 -1)
//    m_max: maximum order of the spherical harmonics expansion (0...n_max). m_max = 0 for the zonal harmonics only
//
// Return value:
//    acceleration vector in the inertial reference frame [m/s^2]
// ----------------------------------------------
Vector AccelHarmonic(const Vector& r, const Matrix& E, int n_max, int m_max);

// ----------------------------------------------
// gravField
// ----------------------------------------------
// Purpose:
//    Computes the differential equation for the gravitational field
//
// Parameters:
//    n: number of equations
//    t: initial time of the integration
//    x: initial state vector
//    f: output vector
//    param: pointer to the structure containing the parameters of the gravitational field
//
// Return value:
//    0: success
//    1: failure
// ----------------------------------------------
int gravField(int n, double t, double x[], double f[], void* param);

// ----------------------------------------------
// integrateOrbit
// ----------------------------------------------
// Purpose:
//    Integrates the differential equation for the gravitational field and computes the position and velocity vectors as well as the orbital elements of the 'new satellite' at the end of the integration
//
// Parameters:
//    s0: initial and final satellite
//    mjd_tt_0: initial time of the integration
//    T: integration time [s]
//    h: step size of integration
//    tol: tolerance for the integration
//    maxNumSteps: maximum number of steps of integration
//    param: pointer to the structure containing the parameters of the gravitational field
//
// Return value:
//    0: success
//    1: failure
// ----------------------------------------------
int integrateOrbit(Satellite& s0, double T, double h, double tol, int maxNumSteps, void* param);
#endif  // FORCE_H