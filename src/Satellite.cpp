#include "../include/Satellite.h"

#include <cmath>
#include <cstdarg>  // in order to use an arbitrary number of arguments
#include <iomanip>  // in order to use std::setprecision
#include <iostream>
#include <sstream>    // in order to use std::istringstream
#include <stdexcept>  // in order to use exceptions

#include "../include/LinearAlgebra.h"
#include "../include/ReferenceSystems.h"
#include "../include/constants.h"
#include "../include/kepler.h"
#include "../include/newton.h"

Satellite::Satellite(const std::string& line0, const std::string& line1, const std::string& line2) : sat_name(line0), line1(line1), line2(line2) {
  // Read TLE data
  std::istringstream iss1(line1);
  std::istringstream iss2(line2);
  std::string temp;
  iss1 >> temp >> temp >> temp >> temp >> dn;
  iss2 >> temp >> temp >> i >> Omega >> temp >> omega >> M;

  // Convert TLE data to useful values
  // Line 1
  std::istringstream iss_year(line1.substr(18, 2));
  iss_year >> epoch_year;
  if (epoch_year < 57)  // 1957 is the first year of the GPS system
    epoch_year += 2000;
  else
    epoch_year += 1900;
  std::istringstream(line1.substr(20, 12)) >> epoch_day;
  std::istringstream(line1.substr(44, 1) + "0." + line1.substr(45, 5) + "e" + line1.substr(50, 2)) >> ddn;  // scientific notation
  // Line 2
  std::istringstream("0." + line2.substr(26, 7)) >> e;
  std::istringstream(line2.substr(52, 11)) >> n;
  std::istringstream(line2.substr(63, 5)) >> rev_num;

  // Convert data to SI units
  i *= DEG_TO_RAD;
  Omega *= DEG_TO_RAD;
  omega *= DEG_TO_RAD;
  M *= DEG_TO_RAD;
  dn *= REV_TO_RAD / (DAY_TO_SEC * DAY_TO_SEC);
  ddn *= REV_TO_RAD / (DAY_TO_SEC * DAY_TO_SEC * DAY_TO_SEC);
  n *= REV_TO_RAD / DAY_TO_SEC;

  a = pow(GM_EARTH / (n * n), 1.0 / 3.0);  // semimajor axis

  solve_kepler_equation();                             // Solve kepler equation to get the eccentric anomaly
  nu = atan(sqrt(1 - e * e) * sin(E) / (cos(E) - e));  // True anomaly
  nu += 2 * PI * ((nu < 0) ? 1 : 0);                   // Make sure that nu is in the range [0, 2PI]
  set_position_velocity();                             // Set position and velocity vectors
  get_gregorian_date(epoch_year, epoch_day, MM, D);
  mjd_TT = get_mjd_TT(epoch_year, MM, D);
  r_GCRF = GCRF2Perifocal(i, Omega, omega) * r_peri;
  v_GCRF = GCRF2Perifocal(i, Omega, omega) * v_peri;
}

void Satellite::print() const {
  std::cout << sat_name << std::endl
            << std::endl;
  std::cout << std::setprecision(12);
  // std::cout << "Epoch year\t\t\t\t\t" << epoch_year << std::endl;
  // std::cout << "Epoch day\t\t\t\t\t" << epoch_day << std::endl;
  // std::cout << "Epoch month\t\t\t\t\t" << MM << std::endl;
  // std::cout << "Epoch day of month\t\t\t\t" << D << std::endl;
  // std::cout << "Epoch MJD TT\t\t\t\t\t" << mjd_TT << std::endl;
  // std::cout << "First time derivative of the mean motion\t" << dn << std::endl;
  // std::cout << "Second time derivative of the mean motion\t" << ddn << std::endl;
  std::cout << "Inclination\t\t\t\t\t" << i << std::endl;
  std::cout << "Right ascension (Omega)\t\t\t\t" << Omega << std::endl;
  std::cout << "Eccentricity\t\t\t\t\t" << e << std::endl;
  std::cout << "Argument of perigee (omega)\t\t\t" << omega << std::endl;
  std::cout << "Mean anomaly\t\t\t\t\t" << M << std::endl;
  std::cout << "True anomaly\t\t\t\t\t" << nu << std::endl;
  std::cout << "Eccentric anomaly\t\t\t\t" << E << std::endl;
  std::cout << "Mean motion\t\t\t\t\t" << n << std::endl;
  std::cout << "Revolution number at epoch\t\t\t" << rev_num << std::endl;
  std::cout << "Semimajor axis\t\t\t\t\t" << a << std::endl;
}

// Vector Satellite::get_position_GCRF() const {
//   // Rotation matrix from perifocal to GCRF
//   Matrix R = GCRF2Perifocal(i, Omega, omega);
//   return R.transpose() * r_peri;
// }

// Vector Satellite::get_position_ITRF() const {
//   // Rotation matrix from GCRF to ITRF
//   Matrix R = GCRF2ITRF(mjd_TT);
//   return R * get_position_GCRF();
// }

int Satellite::solve_kepler_equation() {
  // Initial guess (see Montenbruck & Gill, 2000, p. 24)
  if (e > 0.8)
    E = PI;
  else
    E = M;
  double tol = 1e-12;
  int max_iter = 100;
  args_kepler args = {e, M};

  Vector v;
  // Solve kepler equation
  if (newton(&E, tol, max_iter, kepler_equation, diff_kepler_equation, &args)) {
    throw std::runtime_error("Newton method failed to converge");
    return 1;
  }
  // std::cout << "f = " << kepler_equation(3, M, e) << ", f' = " << diff_kepler_equation(3, M, e) << std::endl;
  return 0;
}

void Satellite::set_position_velocity() {
  double x = a * (cos(E) - e);
  double y = a * sqrt(1 - e * e) * sin(E);
  double z = 0;

  double vx = -sqrt(GM_EARTH / a) * sin(E) / (1 - e * cos(E));
  double vy = sqrt(GM_EARTH / a) * sqrt(1 - e * e) * cos(E) / (1 - e * cos(E));
  double vz = 0;

  r_peri = Vector(x, y, z);
  v_peri = Vector(vx, vy, vz);
}

void Satellite::set_orbital_elements(Vector r, Vector v) {
  Vector h = r.cross(v);
  Vector normal = Vector(0, 0, 1).cross(h);
  Vector e_vec = (r * (v.norm() * v.norm() - GM_EARTH / r.norm()) - v * r.dot(v)) / GM_EARTH;
  double p = h.norm() * h.norm() / GM_EARTH;
  Vector W = h.normalized();
  i = atan(sqrt(W(0) * W(0) + W(1) * W(1)) / W(2));
  printf("W = (%lf, %lf, %lf)\n", W(0), W(1), W(2));
  Omega = atan(-W(0) / W(1));
  double Omega1 = acos(normal(0) / normal.norm());
  double omega1 = acos(normal.dot(e_vec) / (normal.norm() * e_vec.norm()));
  std::cout << "Omega1 = " << Omega1 << std::endl;
  a = 1. / (2. / r.norm() - v.norm() * v.norm() / GM_EARTH);
  std::cout << "omega1 = " << omega1 << std::endl;
  n = sqrt(GM_EARTH / (a * a * a));
  e = sqrt(1 - p / a);
  E = atan(r.dot(v) / (a * a * n) / (1. - r.norm() / a));
  M = E - e * sin(E);

  double u = atan(r(2) / (-r(0) * W(1) + r(1) * W(0)));
  nu = atan(sqrt(1 - e * e) * sin(E) / (cos(E) - e));  // True anomaly
  omega = u - nu;

  // Make sure all angles are in [0, 2PI]
  i += 2 * PI * ((i < 0) ? 1 : 0);
  Omega += 2 * PI * ((Omega < 0) ? 1 : 0);
  E += 2 * PI * ((E < 0) ? 1 : 0);
  M += 2 * PI * ((M < 0) ? 1 : 0);
  nu += 2 * PI * ((nu < 0) ? 1 : 0);
  omega += 2 * PI * ((omega < 0) ? 1 : 0);
}