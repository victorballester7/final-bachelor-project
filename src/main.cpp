// Finally, here is the `main.cpp` file that creates a `TLE` object and prints out the satellite's position and velocity at the epoch time:

#include <stdlib.h>

#include <fstream>
#include <iostream>

#include "../include/LinearAlgebra.h"
#include "../include/Satellite.h"
#include "../include/force.h"
#include "../include/misc.h"

int main(int argc, char const *argv[]) {
  // read TLE data from data/tle/stations.txt
  if (argc != 2) {
    std::cout << "Usage: ./main <satellite>" << std::endl;
    return 1;
  }
  std::string line0, line1, line2;
  std::string satellite = argv[1];
  std::string filename = "data/tle/" + satellite + ".txt";
  std::ifstream file_in(filename);

  if (!file_in.is_open()) {
    std::cout << "Error opening file" << std::endl;
    return 1;
  }
  // int LINE = 0;
  // // we are not interested in the first 66 lines
  // for (int i = 0; i < LINE - 1; i++)
  //   std::getline(file_in, line0);

  // std::getline(file_in, line1);
  // std::getline(file_in, line2);
  // Satellite tle("NUTSAT", line1, line2);
  // tle.print();

  // tle.set_orbital_elements(tle.r_ECI, tle.v_ECI);

  // tle.print();
  // This is for checking the data of the orbital elements for all the satellites in the file_in
  // ---------------------------------------------------------------
  // int n_lines = numLines(filename);
  // printf("Number of lines: %d\n", n_lines);
  // for (int j = 0; j < n_lines; j += 2) {
  //   // std::getline(file_in, line0);
  //   std::getline(file_in, line1);
  //   std::getline(file_in, line2);

  //   Satellite tle("NUTSAT" + std::to_string(j / 2), line1, line2);
  //   std::string name = tle.sat_name;
  //   double e, a, i, Omega, omega, M, n, nu, E;
  //   e = tle.e;
  //   a = tle.a;
  //   i = tle.i;
  //   Omega = tle.Omega;
  //   omega = tle.omega;
  //   M = tle.M;
  //   n = tle.n;
  //   nu = tle.nu;
  //   E = tle.E;

  //   tle.set_orbital_elements(tle.r_ECI, tle.v_ECI);

  //   // std::cout << "Satellite name: " << name << " (line " << (j + 1) << ")" << std::endl;
  //   std::cout << "Satellite name: " << name << std::endl;
  //   std::cout << "Line: " << (j + 1) << std::endl;
  //   std::cout << "Eccentricity: \t\t\t\t" << (fabs(e - tle.e) < 1.e-8) << std::endl;
  //   std::cout << "Semimajor axis: \t\t\t" << (fabs(a - tle.a) < 1.e-8) << std::endl;
  //   std::cout << "Inclination: \t\t\t\t" << (fabs(i - tle.i) < 1.e-8) << std::endl;
  //   std::cout << "Right ascension of the ascending node:\t" << (fabs(Omega - tle.Omega) < 1.e-8) << std::endl;
  //   std::cout << "Argument of perigee:\t\t\t" << (fabs(omega - tle.omega) < 1.e-8) << std::endl;
  //   std::cout << "Mean anomaly: \t\t\t\t" << (fabs(M - tle.M) < 1.e-8) << std::endl;
  //   std::cout << "Mean motion: \t\t\t\t" << (fabs(n - tle.n) < 1.e-8) << std::endl;
  //   std::cout << "True anomaly: \t\t\t\t" << (fabs(nu - tle.nu) < 1.e-8) << std::endl;
  //   std::cout << "Eccentric anomaly: \t\t\t" << (fabs(E - tle.E) < 1.e-8) << std::endl;
  //   std::cout << std::endl;
  //   // tle.print();
  // }
  // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -

  // /*
  int n_lines = numLines(filename);
  std::ofstream file_out, file_out2;
  file_out.open("data/tle/errors_" + satellite + ".txt");
  file_out2.open("data/tle/orbit_" + satellite + ".txt");
  if (!file_out.is_open() || !file_out2.is_open()) {
    std::cout << "Error opening file" << std::endl;
    return 1;
  }
  std::getline(file_in, line1);
  std::getline(file_in, line2);
  Satellite s(satellite, line1, line2);
  std::cout << "Position satellite " << 0 << ": " << s.r_ECI << std::endl;

  file_in.clear();
  file_in.seekg(0, std::ios::beg);
  // printf("\n\n");
  args_gravField args = {8, 8, s.mjd_TT, false, false, false};  // n_max, m_max, pointEarth, sun, moon
  double tol = 1e-8;
  double h, T, t0 = s.mjd_TT, error;
  int maxNumSteps = 10000;
  int steps = 100;
  int numTLEs = 15;
  // comparing Omega
  // for (int i = 0; i < n_lines / 2; i++) {
  for (int i = 0; i < numTLEs; i++) {
    std::getline(file_in, line1);
    std::getline(file_in, line2);

    Satellite tle(satellite + std::to_string(i), line1, line2);
    T = (tle.mjd_TT - s.mjd_TT) * 86400;
    if (T < 1.e-2) continue;
    h = T / steps;
    if (integrateOrbit(s, T, h, tol, maxNumSteps, &args)) {
      std::cout << "Integration failed" << std::endl;
      return 1;
    }
    // file_out << s.mjd_TT - t0 << " " << tle.Omega << " " << s.Omega << std::endl;
    // file_out << s.i << " " << tle.Omega << " " << s.Omega << std::endl;
    // std::cout << s.r_ECI(0) << " " << s.r_ECI(1) << " " << s.r_ECI(2) << std::endl;
    // std::cout << time_posXYZ(i, 1) << " " << time_posXYZ(i, 2) << " " << time_posXYZ(i, 3) << std::endl;
    // for (int j = 0; j < 3; j++) error += fabs(s.r_ECI(j) - tle.r_ECI(j));                            // error in L1 norm
    for (int j = 0; j < 3; j++) error += (s.r_ECI(j) - tle.r_ECI(j)) * (s.r_ECI(j) - tle.r_ECI(j));  // error in L2 norm
    error = sqrt(error);
    std::cout << "Position satellite " << i << ": " << tle.r_ECI << std::endl;
    file_out << s.mjd_TT - t0 << " " << error << std::endl;
  }
  // */

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