#ifndef FORCE_H
#define FORCE_H

#include "LinearAlgebra.h"
#include "Satellite.h"

typedef struct {
  int n_max;
  int m_max;
} args_gravField;

// ----------------------------------------------
// AccelHarmonic
// ----------------------------------------------
// Purpose:
//    Calculate acceleration due to gravity using the spherical harmonics expansion of the gravitational potential
//
// Parameters:
//    r: position vector in the inertial reference frame [m]
//    E: transformation matrix from the inertial reference frame to the reference frame of the spherical harmonics expansion (fixed to the Earth)
//    GM: gravitational parameter of the central body [m^3/s^2]
//    R_ref: reference radius of the spherical harmonics expansion [m]
//    CS: cosine and sine coefficients of the spherical harmonics expansion
//    n_max: maximum degree of the spherical harmonics expansion (0...N_JDM3 -1)
//    m_max: maximum order of the spherical harmonics expansion (0...n_max). m_max = 0 for the zonal harmonics only
//
// Return value:
//    acceleration vector in the inertial reference frame [m/s^2]
// ----------------------------------------------
Vector AccelHarmonic(const Vector& r, const Matrix& E, double GM, double R_ref, const Matrix& CS, int n_max, int m_max);

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
//    s0: initial satellite
//    sT: final satellite
//    T: integration time
//    hmin: minimum step size of integration
//    hmax: maximum step size of integration
//    tol: tolerance for the integration
//    maxNumSteps: maximum number of steps of integration
//    n_max: maximum degree of the spherical harmonics expansion (0...N_JDM3 -1)
//    m_max: maximum order of the spherical harmonics expansion (0...n_max). m_max = 0 for the zonal harmonics only
//
// Return value:
//    0: success
//    1: failure
// ----------------------------------------------
int integrateOrbit(Satellite& s0, Satellite& sT, double T, double hmin, double hmax, double tol, int maxNumSteps, int n_max, int m_max);

#endif  // FORCE_H