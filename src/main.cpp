// Finally, here is the `main.cpp` file that creates a `TLE` object and prints out the satellite's position and velocity at the epoch time:

#include <fstream>
#include <iostream>

#include "../include/LinearAlgebra.h"
#include "../include/Satellite.h"
#include "../include/force.h"
int main() {
  // read TLE data from /data/stations.txt
  std::string line0, line1, line2;
  std::ifstream file("data/stations.txt");

  if (!file.is_open()) {
    std::cout << "Error opening file" << std::endl;
    return 1;
  }
  int LINE = 67;
  // we are not interested in the first 66 lines
  for (int i = 0; i < LINE - 1; i++)
    std::getline(file, line0);

  std::getline(file, line0);
  std::getline(file, line1);
  std::getline(file, line2);
  Satellite tle(line0, line1, line2);
  tle.print();

  tle.set_orbital_elements(tle.r_GCRF, tle.v_GCRF);

  tle.print();

  // args_gravField numHarmonics = {8, 8};  // n_max, m_max
  // Satellite s_forw, s_back;
  // double h_min = 0.01;
  // double h_max = 0.05;
  // double tol = 1e-12;
  // int maxNumSteps = 100;
  // integrateOrbit(tle, s_forw, 1, h_min, h_max, tol, maxNumSteps, numHarmonics.n_max, numHarmonics.m_max);

  // integrateOrbit(tle, s_back, -1, h_min, h_max, tol, maxNumSteps, numHarmonics.n_max, numHarmonics.m_max);

  // printf("\n\n");
  // s_forw.print();

  // printf("\n\n");
  // s_back.print();

  // std::cout << "Mean motion (forward): " << s_forw.M << std::endl;
  // std::cout << "Mean motion (backward): " << s_back.M << std::endl;

  // double diff1M = (s_forw.M - s_back.M) / (2 * tol);
  // double diff2M = (s_forw.M - 2 * tle.M + s_back.M) / (tol * tol);

  // std::cout << "First derivative mean motion: " << diff1M << std::endl;
  // std::cout << "Second derivative mean motion: " << diff2M << std::endl;

  return 0;
}