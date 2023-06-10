// Finally, here is the `main.cpp` file that creates a `TLE` object and prints out the satellite's position and velocity at the epoch time:

#include <stdlib.h>

#include <fstream>
#include <iostream>
#include <sstream>  // in order to use std::istringstream

#include "../include/LinearAlgebra.h"
#include "../include/ReferenceSystems.h"
#include "../include/Satellite.h"
#include "../include/force.h"
#include "../include/misc.h"

int main(int argc, char const *argv[]) {
  // read TLE data from data/tle/stations.txt
  if (argc != 2) {
    std::cout << "Usage: ./main <satellite>" << std::endl;
    return 1;
  }
  std::string satellite = argv[1];
  // Read TEME data from file
  // std::string filename_in = "data/teme/" + satellite + "_tles.txt";
  // std::string filename_in = "data/teme/" + satellite + "_1h_interval.txt";
  std::string filename_in = "data/teme/" + satellite + "_1min_interval.txt";
  std::string filename_out = "data/errors/" + satellite + "_errors.txt";
  std::string filename_out_orbit = "data/orbit/" + satellite + "_orbit.txt";
  std::ifstream file_in(filename_in);
  std::ofstream file_out(filename_out), file_out_orbit(filename_out_orbit);
  if (!file_in.is_open() || !file_out.is_open() || !file_out_orbit.is_open()) {
    std::cout << "Error opening file" << std::endl;
    return 1;
  }

  int TLEs = 1000;
  std::string line;
  double mjd_utc;
  double x, y, z, vx, vy, vz;
  std::getline(file_in, line);
  std::istringstream iss(line);
  iss >> mjd_utc >> x >> y >> z >> vx >> vy >> vz;
  // Reset file pointer to beginning of file
  file_in.clear();
  file_in.seekg(0, std::ios::beg);
  Vector r_teme = Vector(x, y, z);
  Vector v_teme = Vector(vx, vy, vz);
  Satellite s = Satellite(mjd_utc, r_teme, v_teme);
  r_teme.print();
  v_teme.print();
  s.r_ECI.print();
  s.v_ECI.print();

  args_gravField args = {8, 8, s.mjd_TT, true, false, false, false, false, 10000., 450000.};  // n_max, m_max, initial mjd_TT, pointEarth, sun, moon, solarRad, atmo_drag, area, mass
  double tol = 1e-12;
  double h, T, t0 = s.mjd_TT, error, temp, t1;
  int maxNumSteps = 1000;
  int inter_steps = 100;
  for (int i = 0; i < TLEs; i++) {
    std::getline(file_in, line);
    std::istringstream iss(line);
    iss >> mjd_utc >> x >> y >> z >> vx >> vy >> vz;
    r_teme = Vector(x, y, z);
    v_teme = Vector(vx, vy, vz);
    Satellite t = Satellite(mjd_utc, r_teme, v_teme);
    // set precision
    std::cout.precision(16);
    T = (t.mjd_TT - s.mjd_TT) * 86400;
    std::cout << "T = " << T << std::endl;
    if (T > 1) {  // t0 and t1 are different
      h = T / inter_steps;
      if (integrateOrbit(s, T, h, tol, maxNumSteps, &args)) {
        std::cout << "Integration failed" << std::endl;
        return 1;
      }
    }
    file_out_orbit << s.r_ECI(0) << " " << s.r_ECI(1) << " " << s.r_ECI(2) << " " << s.v_ECI(0) << " " << s.v_ECI(1) << " " << s.v_ECI(2) << " " << s.r_ECI.norm() - R_JGM3 << std::endl;

    error = (s.r_ECI - t.r_ECI).norm();
    file_out << t.mjd_TT - t0 << " " << error << std::endl;

    // // plot only error in z
    // file_out << t1 - t0 << " " << s.r_ECI.norm() - R_JGM3 << std::endl;
  }
  // int LINE = 0;
  // // we are not interested in the first 66 lines
  // for (int i = 0; i < LINE - 1; i++)
  //   std::getline(file_in, line0);

  // Satellite s = Satellite("prova", "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  47533", "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.8241915741366");

  // double time = 0. * 86400.;
  // int steps = 1000;
  // double h = time / steps;
  // double tol = 1e-8;
  // int maxNumSteps = 10000;
  // args_gravField args = {8, 8, s.mjd_TT, false, false, false};  // n_max, m_max, pointEarth, sun, moon

  // if (integrateOrbit(s, time, h, tol, maxNumSteps, &args)) {
  //   std::cout << "Integration failed" << std::endl;
  //   return 1;
  // }

  // std::cout << "mjd_TT: " << s.mjd_TT << std::endl;

  // std::cout << "r_ECI: " << s.r_ECI << std::endl;
  // std::cout << "v_ECI: " << s.v_ECI << std::endl;
  // std::cout << "r_TEME: " << ECI2TEME(s.mjd_TT) * s.r_ECI << std::endl;
  // std::cout << "v_TEME: " << ECI2TEME(s.mjd_TT) * s.v_ECI << std::endl;

  // int n_lines = numLines(filename);
  // // comparing Omega
  // // for (int i = 0; i < n_lines / 2; i++) {
  // for (int i = 0; i < n_lines / 2; i++) {
  //   std::getline(file_in, line1);
  //   std::getline(file_in, line2);

  //   Satellite tle(satellite + std::to_string(i), line1, line2);
  //   std::cout << "Position satellite " << i << ": " << tle.r_ECI << std::endl;
  // }

  // for (int j = 0; j < 3; j++) error += fabs(s.r_ECI(j) - time_posXYZ(i, j + 1));
  // file_out << s.mjd_TT - t0 << " " << error << std::endl;
  // while (s.mjd_TT < time_posXYZ(numTLEs - 1, 0) - tol) {
  //   i++;
  //   error = 0;
  //   T = (time_posXYZ(i, 0) - s.mjd_TT) * 86400;
  //   h = T / steps;
  //   if (integrateOrbit(s, T, h, tol, maxNumSteps, &args)) {
  //     std::cout << "Integration failed" << std::endl;
  //     return 1;
  //   }
  //   // std::cout << s.mjd_TT << " " << time_posXYZ(i, 0) << std::endl;

  //   // printf("%d %f %f\n", i, s.mjd_TT, error);
  // }

  return 0;
}