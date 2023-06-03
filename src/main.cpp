// Finally, here is the `main.cpp` file that creates a `TLE` object and prints out the satellite's position and velocity at the epoch time:

#include <stdlib.h>

#include <fstream>
#include <iostream>

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
  //
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // ATEEEEEEEEEEEEEEEEEEEENCCCCCCCCCCCCCCCCCCIIIIIIIIIIIOOOOOOOOOOOOOO QUE S'HA DE CANVIAR EL mdj_tt per el mjd_ut1 en les funcions GMST dins the GAST i RotationMatrix.
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // April 6, 2004, 7:51:27.946047 UT1
  // double mjd_ut1 = get_mjd(2004, 4, 6.3274067829);
  // double mjd_tt = get_mjd(2004, 4, 6.328154745);
  // std::cout << "MJD: " << mjd_ut1 << std::endl;
  // double gmst = GMST(mjd_ut1);
  // double gast = GAST(mjd_tt);

  // std::cout << "GMST: " << gmst / DEG_TO_RAD << std::endl;
  // std::cout << "GAST: " << gast / DEG_TO_RAD << std::endl;
  // // double mjd_tt = 53101.328154745;
  // // Vector r = Vector(-1033.4793830, 7901.2952754, 6380.3565958);
  // // Vector v = Vector(-3.225636520, -2.872451450, 5.531924446);
  // Vector r = Vector(-1033.4750313, 7901.3055856, 6380.3445328);

  // double dpsi = 0.0;
  // double deps = 0.0;

  // NutationAngles(mjd_tt, dpsi, deps);

  // std::cout << "dpsi: " << dpsi / DEG_TO_RAD << std::endl;
  // std::cout << "deps: " << deps / DEG_TO_RAD << std::endl;

  // Matrix T = J20002ECEF(mjd_tt);

  // Vector r_J2000 = T.transpose() * r;
  // // Vector v_J2000 = T.transpose() * v;

  // Vector r_rot = RotationMatrix(mjd_tt).transpose() * r;

  // std::cout << "r_rot: " << r_rot << std::endl;

  // Vector r_tod_real = Vector(5094.5147804, 6127.3664612, 6380.3445328);
  // Vector v_tod_real = Vector(-4.746088567, 0.786077222, 5.531931288);

  // Vector r_mod_real = Vector(5094.0283745, 6127.8708164, 6380.2485164);
  // Vector v_mod_real = Vector(-4.746263052, 0.786014045, 5.531790562);

  // Vector r_mod = NutationMatrix(mjd_tt).transpose() * r_tod_real;
  // Vector v_mod = NutationMatrix(mjd_tt).transpose() * v_tod_real;

  // // Vector r_gcrf = PrecessionMatrix(mjd_tt).transpose() * r_mod_real;
  // // Vector v_gcrf = PrecessionMatrix(mjd_tt).transpose() * v_mod_real;

  // Vector r_gcrf = PrecessionMatrix(mjd_tt).transpose() * NutationMatrix(mjd_tt).transpose() * r_tod_real;
  // Vector v_gcrf = PrecessionMatrix(mjd_tt).transpose() * NutationMatrix(mjd_tt).transpose() * v_tod_real;

  // std::cout << "r_mod: " << r_mod << std::endl;
  // std::cout << "v_mod: " << v_mod << std::endl;

  // std::cout << "r_gcrf: " << r_gcrf << std::endl;
  // std::cout << "v_gcrf: " << v_gcrf << std::endl;

  // std::cout << "r_J2000: " << r_J2000 << std::endl;

  Satellite s = Satellite("prova", "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  47533", "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.8241915741366");

  double time = 0. * 86400.;
  int steps = 1000;
  double h = time / steps;
  double tol = 1e-8;
  int maxNumSteps = 10000;
  args_gravField args = {8, 8, s.mjd_TT, false, false, false};  // n_max, m_max, pointEarth, sun, moon

  if (integrateOrbit(s, time, h, tol, maxNumSteps, &args)) {
    std::cout << "Integration failed" << std::endl;
    return 1;
  }

  std::cout << "mjd_TT: " << s.mjd_TT << std::endl;

  std::cout << "r_ECI: " << s.r_ECI << std::endl;
  std::cout << "v_ECI: " << s.v_ECI << std::endl;
  std::cout << "r_TEME: " << ECI2TEME(s.mjd_TT) * s.r_ECI << std::endl;
  std::cout << "v_TEME: " << ECI2TEME(s.mjd_TT) * s.v_ECI << std::endl;

  // /*
  // int n_lines = numLines(filename);
  // std::ofstream file_out, file_out2;
  // file_out.open("data/tle/errors_" + satellite + ".txt");
  // file_out2.open("data/tle/orbit_" + satellite + ".txt");
  // if (!file_out.is_open() || !file_out2.is_open()) {
  //   std::cout << "Error opening file" << std::endl;
  //   return 1;
  // }
  // std::getline(file_in, line1);
  // std::getline(file_in, line2);
  // Satellite s(satellite, line1, line2);
  // std::cout << "Position satellite " << 0 << ": " << s.r_ECI << std::endl;

  // file_in.clear();
  // file_in.seekg(0, std::ios::beg);
  // // printf("\n\n");
  // args_gravField args = {8, 8, s.mjd_TT, false, false, false};  // n_max, m_max, pointEarth, sun, moon
  // double tol = 1e-8;
  // double h, T, t0 = s.mjd_TT, error;
  // int maxNumSteps = 10000;
  // int steps = 100;
  // int numTLEs = 15;
  // // comparing Omega
  // // for (int i = 0; i < n_lines / 2; i++) {
  // for (int i = 0; i < numTLEs; i++) {
  //   std::getline(file_in, line1);
  //   std::getline(file_in, line2);

  //   Satellite tle(satellite + std::to_string(i), line1, line2);
  //   T = (tle.mjd_TT - s.mjd_TT) * 86400;
  //   if (T < 1.e-2) continue;
  //   h = T / steps;
  //   if (integrateOrbit(s, T, h, tol, maxNumSteps, &args)) {
  //     std::cout << "Integration failed" << std::endl;
  //     return 1;
  //   }
  //   // file_out << s.mjd_TT - t0 << " " << tle.Omega << " " << s.Omega << std::endl;
  //   // file_out << s.i << " " << tle.Omega << " " << s.Omega << std::endl;
  //   // std::cout << s.r_ECI(0) << " " << s.r_ECI(1) << " " << s.r_ECI(2) << std::endl;
  //   // std::cout << time_posXYZ(i, 1) << " " << time_posXYZ(i, 2) << " " << time_posXYZ(i, 3) << std::endl;
  //   // for (int j = 0; j < 3; j++) error += fabs(s.r_ECI(j) - tle.r_ECI(j));                            // error in L1 norm
  //   for (int j = 0; j < 3; j++) error += (s.r_ECI(j) - tle.r_ECI(j)) * (s.r_ECI(j) - tle.r_ECI(j));  // error in L2 norm
  //   error = sqrt(error);
  //   std::cout << "Position satellite " << i << ": " << tle.r_ECI << std::endl;
  //   file_out << s.mjd_TT - t0 << " " << error << std::endl;
  // }
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