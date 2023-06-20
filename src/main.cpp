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
  if (argc != 7) {
    std::cout << "Usage: ./main <satellite> <t/f pointEarth> <t/f sun> <t/f moon> <t/f solarRad> <t/f atmoDrag>" << std::endl;
    return 1;
  }
  std::string satellite = argv[1];
  std::string filename_in_bstar = "data/bstar/" + satellite + ".txt";
  std::ifstream file_in_bstar(filename_in_bstar);
  std::string line;
  double bstar, Am;

  if (!file_in_bstar.is_open()) {
    bstar = 0.;
    Am = 0.01;  // TDRS-3
    // Am = 0.075;
  } else {
    std::getline(file_in_bstar, line);
    std::istringstream iss(line);
    iss >> bstar;
    file_in_bstar.close();
    Am = bstar * 2 / (0.157 * C_D);
  }

  printf("Am = %f\n", Am);

  args_gravField args = {8, 8, false, false, false, false, false, Am};  // n_max, m_max, initial mjd_TT, pointEarth, sun, moon, solarRad, atmo_drag, area_mass

  if (argv[2][0] == 't') args.pointEarth = true;
  if (argv[3][0] == 't') args.sun = true;
  if (argv[4][0] == 't') args.moon = true;
  if (argv[5][0] == 't') args.solar_rad = true;
  if (argv[6][0] == 't') args.atmo_drag = true;

  std::string end_name = "";

  if (args.pointEarth)
    end_name += "_pointEarth";
  else
    end_name += "_sphHarm";

  if (args.sun)
    end_name += "_sun";
  if (args.moon)
    end_name += "_moon";
  if (args.solar_rad)
    end_name += "_solarRad";
  if (args.atmo_drag)
    end_name += "_atmoDrag";

  end_name += ".txt";

  // Read TEME data from file
  // std::string filename_in = "data/teme/" + satellite + "_tles.txt";
  // std::string filename_in = "data/teme/" + satellite + "_1h_interval.txt";
  std::string filename_in = "data/teme/" + satellite + "_1min_interval.txt";
  std::string filename_out = "data/errors/" + satellite + end_name;
  std::string filename_out_orbit = "data/orbit/" + satellite + end_name;
  std::ifstream file_in(filename_in);
  std::ofstream file_out(filename_out), file_out_orbit(filename_out_orbit);
  if (!file_in.is_open() || !file_out.is_open() || !file_out_orbit.is_open()) {
    std::cout << "Error opening file" << std::endl;
    return 1;
  }

  std::cout << "This is a SAAAAAAAAAAAAAAMPLEEEEEEEEEEEEEEE" << std::endl;
  std::cout.precision(16);
  std::cout << Sun(UTC2TT(53827.)) << std::endl;

  int TLEs = 10000;
  double mjd_utc;
  double x, y, z, vx, vy, vz, errTLE;
  std::getline(file_in, line);
  std::istringstream iss(line);
  iss >> mjd_utc >> x >> y >> z >> vx >> vy >> vz;
  // Reset file pointer to beginning of file
  file_in.clear();
  file_in.seekg(0, std::ios::beg);
  Vector r_teme = Vector(x, y, z);
  Vector v_teme = Vector(vx, vy, vz);
  Satellite s = Satellite(mjd_utc, r_teme, v_teme);

  double tol = 1e-12;
  double h, T, t0 = s.mjd_TT, error, temp, t1, height;
  int maxNumSteps = 1000;
  int inter_steps = 1;
  double Errors[TLEs];
  std::cout << "t0 = " << t0 << std::endl;
  for (int i = 0; i < TLEs; i++) {
    std::getline(file_in, line);
    std::istringstream iss(line);
    iss >> mjd_utc >> x >> y >> z >> vx >> vy >> vz >> errTLE;
    r_teme = Vector(x, y, z);
    v_teme = Vector(vx, vy, vz);
    Satellite t = Satellite(mjd_utc, r_teme, v_teme);
    // set precision
    std::cout.precision(16);
    T = (t.mjd_TT - s.mjd_TT) * 86400;
    // std::cout << "i = " << i << " "
    //           << "T = " << T << std::endl;
    if (T > 1) {  // t0 and t1 are different
      h = T / inter_steps;
      if (integrateOrbit(s, T, h, tol, maxNumSteps, &args)) {
        std::cout << "Integration failed" << std::endl;
        return 1;
      }
    }

    Vector r_ECI = s.r_ECI / 1000;  // km
    Vector v_ECI = s.v_ECI / 1000;  // km/s

    // height = Height(NutationMatrix(s.mjd_TT) * PrecessionMatrix(s.mjd_TT) * r_ECI);
    file_out_orbit << s.mjd_TT << r_ECI(0) << " " << r_ECI(1) << " " << r_ECI(2) << " " << v_ECI(0) << " " << v_ECI(1) << " " << v_ECI(2) << " " << r_ECI.norm() - R_Earth / 1000 << std::endl;

    error = (s.r_ECI - t.r_ECI).norm() / 1000;
    file_out << t.mjd_TT - t0 << " " << error << " " << errTLE << std::endl;

    Errors[i] = error;
    // // plot only error in z
    // file_out << t1 - t0 << " " << s.r_ECI.norm() - R_Earth << std::endl;
  }

  std::cout << "Variance of data: " << variance(Errors, TLEs) << std::endl;

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