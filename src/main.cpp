// Finally, here is the `main.cpp` file that creates a `TLE` object and prints out the satellite's position and velocity at the epoch time:

#include <fstream>
#include <iostream>

#include "../include/LinearAlgebra.h"
#include "../include/Satellite.h"
#include "../include/force.h"
#include "../include/misc.h"

int main() {
  // read TLE data from data/tle/stations.txt
  std::string line0, line1, line2;
  std::string satellite = "STARLINK";
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

  int n_lines = numLines(filename);
  // printf("Number of lines: %d\n", n_lines);
  std::ofstream file_out;
  file_out.open("data/tle/orbit_" + satellite + ".txt");
  if (!file_out.is_open()) {
    std::cout << "Error opening file" << std::endl;
    return 1;
  }
  file_out << "# real positions of the satellite in ECI coordinates" << std::endl;
  for (int j = 0; j < n_lines; j += 2) {
    // std::getline(file_in, line0);
    std::getline(file_in, line1);
    std::getline(file_in, line2);

    Satellite tle("NUTSAT" + std::to_string(j / 2), line1, line2);
    file_out << tle.r_ECI(0) << " " << tle.r_ECI(1) << " " << tle.r_ECI(2) << std::endl;
  }
  file_out << "\n\n# integrated positions of the satellite in ECI coordinates" << std::endl;

  args_gravField args = {8, 8, true, false, false};  // n_max, m_max, pointEarth, sun, moon
  double tol = 1e-8;
  double h = 100;
  int maxNumSteps = 10000;
  file_in.clear();
  file_in.seekg(0, std::ios::beg);

  std::getline(file_in, line1);
  std::getline(file_in, line2);
  Satellite s("NUTSAT", line1, line2);
  s.print();
  file_out << s.r_ECI(0) << " " << s.r_ECI(1) << " " << s.r_ECI(2) << std::endl;
  for (int j = 2; j < n_lines; j += 2) {
    // std::getline(file_in, line0);
    std::getline(file_in, line1);
    std::getline(file_in, line2);

    Satellite tmp("NUTSAT", line1, line2);

    // std::cout << "time of integration: " << tmp.mjd_TT - s.mjd_TT << std::endl;

    if (integrateOrbit(s, (tmp.mjd_TT - s.mjd_TT) * 86400, h, tol, maxNumSteps, &args)) {
      std::cout << "Integration failed" << std::endl;
      return 1;
    }
    file_out << s.r_ECI(0) << " " << s.r_ECI(1) << " " << s.r_ECI(2) << std::endl;
  }

  // args_gravField args = {8, 8, true, false, false};  // n_max, m_max, pointEarth, sun, moon
  // Satellite s_forward, s_backward;
  // double tol = 1e-8;
  // int maxNumSteps = 1000;
  // double time_of_integration = 16000;

  // if (integrateOrbit(tle, s_forward, time_of_integration, tol, maxNumSteps, &args)) {
  //   printf("Integration failed\n");
  //   return 1;
  // }

  // if (integrateOrbit(tle, s_backward, -time_of_integration, tol, maxNumSteps, &args)) {
  //   printf("Integration failed\n");
  //   return 1;
  // }

  // printf("\n\n");
  // s_forward.print();

  // printf("\n\n");
  // s_backward.print();

  // std::cout << "Mean motion (forward): " << s_forward.n << std::endl;
  // std::cout << "Mean motion (backward): " << s_backward.n << std::endl;

  // double diff1M = (s_forward.n - s_backward.n) / (2 * time_of_integration);
  // double diff2M = (s_forward.n - 2 * tle.n + s_backward.n) / (time_of_integration * time_of_integration);

  // std::cout << "First derivative mean motion (integrated): " << diff1M << std::endl;
  // std::cout << "First derivative mean motion ('analytical'): " << tle.dn << std::endl;
  // std::cout << "Second derivative mean motion (integrated): " << diff2M << std::endl;
  // std::cout << "Second derivative mean motion ('analytical'): " << tle.ddn << std::endl;

  return 0;
}