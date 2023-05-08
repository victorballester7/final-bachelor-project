#ifndef SATELLITE_H
#define SATELLITE_H

#include <string>

#include "LinearAlgebra.h"
#include "misc.h"
class Satellite {
 public:
  // Constructor
  // line0: satellite name
  // line1: first line of TLE data
  // line2: second line of TLE data
  Satellite() = default;
  Satellite(const std::string& line0, const std::string& line1, const std::string& line2);

  std::string sat_name;
  // Orbital elements
  double i;      // inclination of the orbit [rad]
  double Omega;  // right ascension of the ascending node [rad]
  double e;      // eccentricity of the orbit [rad]
  double omega;  // argument of perigee [rad]
  double M;      // mean anomaly [rad]
  double n;      // mean motion [rad/s]

  double dn;       // first time derivative of the mean motion [rad/s^2]
  double ddn;      // second time derivative of the mean motion [rad/s^3]
  double rev_num;  // revolution number at epoch [rev]

  double nu;  // true anomaly [rad]
  double E;   // eccentric anomaly [rad]
  double a;   // semimajor axis [m]
  double P;   // orbital period [s]

  // Time
  int epoch_year;    // year of the epoch (YYYY)
  double epoch_day;  // day of the epoch in UTC (fractional and starting from 000.DDDDDD) (DDD.DDDDDD)
  int MM;            // month of the epoch (MM)
  double D;          // day of the epoch (DD.DDDDDD)
  double mjd_TT;     // Modified Julian date TT (terrestrial time)
  double epoch_jd;   // Julian date of the epoch (JD)

  // Position and velocity
  Vector r_peri;  // position vector in perifocal coordinates [m]
  Vector v_peri;  // velocity vector in perifocal coordinates [m/s]
  Vector r_ECI;   // position vector in ECI coordinates [m]
  Vector v_ECI;   // velocity vector in ECI coordinates [m/s]

  // Functions
  void print() const;
  void set_orbital_elements(Vector r, Vector v);  // Set orbital elements from position and velocity vectors

 private:
  std::string line1;
  std::string line2;
  int solve_kepler_equation();   // Solve Kepler's equation using Newton's method
  void set_position_velocity();  // Calculate position and velocity vectors in perifocal coordinates
  // Vector get_position_ITRF() const;
  // Vector get_position_GCRF() const;
};

#endif  // SATELLITE_H