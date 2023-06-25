#include "../include/force.h"

#include <cmath>

#include "../include/LinearAlgebra.h"
#include "../include/ReferenceSystems.h"
#include "../include/Satellite.h"
#include "../include/constants.h"
#include "../include/rk78.h"

Vector AccelPointEarth(const Vector& r) {
  double r3 = pow(r.norm(), 3.);
  return r * (-GM_EARTH / r3);
}

Vector AccelPointMass(const Vector& r, const Vector& s, double GM) {
  // relative distance between big body and small body
  Vector r_relSR = s - r;
  double r1 = r_relSR.norm();
  double r2 = r1 * r1;
  double r3 = r2 * r1;
  double s1 = s.norm();
  double s2 = s1 * s1;
  double s3 = s2 * s1;
  // Vector sat_pointMass = r_relSR / r3;
  // Vector earth_pointMass = s / s3;

  // return (r_relSR / r3 - s / s3) * GM;

  // Solution to improve cancellation error (from Vallado: p. 575)
  double Q = (r.dot(r) + 2 * r.dot(r_relSR)) * (s2 + s1 * r1 + r2) / (s3 * r3 * (s1 + r1));
  return (r_relSR * Q - r / s3) * (GM);

  // double q = (r.dot(r) + 2 * r.dot(s)) / s2;
  // double f = q * (3 + q * (3 + q)) / (1 + pow(1 + q, 1.5));
  // return (r - s * f) * (-GM / (s3 * pow(1 + q, 1.5)));

  // static int count = 0;
  // if (count++ < 10) {
  //   std::cout << "value Vallado: " << (sat_pointMass - earth_pointMass) * (GM) << std::endl;
  //   std::cout << "value Roy: " << (r_relSR * Q - r / s3) * GM << std::endl;
  //   count++;
  // }
  // return (sat_pointMass - earth_pointMass) * (GM);
}

Vector AccelSolarRad(const Vector& r, const Vector& r_Sun, double Area_mass) {
  Vector d(3);

  // Relative position vector of spacecraft w.r.t. Sun

  d = r - r_Sun;

  double d3 = pow(d.norm(), 3);  // Distance cubed

  // Acceleration

  // return d * (C_R * Area_mass * P0 * AU * AU / d3);
  return d * (C_R * Area_mass * P0 / d.norm());
}

double Illumination(const Vector& r, const Vector& r_Sun) {
  // From Vallado p. 299
  double alpha_umb, alpha_pen;  // Umbra and penumbra angles
  double zeta, sat_horiz, sat_vert, pen_vert, umb_vert;
  alpha_umb = asin((R_Sun - R_Earth) / r_Sun.norm());
  alpha_pen = asin((R_Sun + R_Earth) / r_Sun.norm());
  if (r.dot(r_Sun) > 0)  // Sunlight
    return 1.;
  zeta = acos(r.dot(r_Sun) / (r.norm() * r_Sun.norm()));
  sat_horiz = r.norm() * cos(zeta);
  sat_vert = r.norm() * sin(zeta);
  pen_vert = R_Earth + sat_horiz * tan(alpha_pen);
  umb_vert = R_Earth - sat_horiz * tan(alpha_umb);
  if (sat_vert < umb_vert)  // Umbra
    return 0.;
  else if (sat_vert > pen_vert)  // Sunlight
    return 1.;
  else  // Partial illumination = Penumbra
    return (sat_vert - umb_vert) / (pen_vert - umb_vert);

  // Montenbruck
  // Vector e_Sun = r_Sun / r_Sun.norm();  // Sun direction unit vector
  // double s = r.dot(e_Sun);              // Projection of s/c position

  // return ((s > 0 || (r - e_Sun * s).norm() > R_Earth) ? 1.0 : 0.0);
}

Vector AccelDrag(const Vector& r, const Vector& v, double mjd_TT, const Matrix& NP, double Area_mass) {
  // Constants

  // Earth angular velocity vector [rad/s]
  const Vector omega = Vector(0.0, 0.0, 7.29212e-5);

  // Variables
  double v_abs, dens;
  Vector v_rel(3), a_tod(3);
  Matrix NP_trans = NP.transpose();

  // Position and velocity in true-of-date system
  Vector r_tod = NP * r;
  Vector v_tod = NP * v;

  // Velocity relative to the Earth's atmosphere
  v_rel = v_tod - omega.cross(r_tod);
  v_abs = v_rel.norm();

  // Atmospheric density due to modified Harris-Priester model
  dens = Density_HP(mjd_TT, r_tod);

  // Acceleration
  a_tod = v_rel * (-0.5 * C_D * Area_mass * dens * v_abs);

  // std::cout << "v_abs^2 = " << v_abs * v_abs << std::endl;
  // std::cout << "dens = " << dens << std::endl;
  // std::cout << "C_D = " << C_D << std::endl;
  // std::cout << "0.5 * C_D * (Area / mass) = " << 0.5 * C_D * (Area / mass) << std::endl;
  // std::cout << "a_tod = " << a_tod.norm() << std::endl;

  return NP_trans * a_tod;
}

double Density_HP(double mjd_TT, const Vector& r_NP) {
  // Constants

  const double upper_limit = 1000.0;  // Upper height limit [km]
  const double lower_limit = 100.0;   // Lower height limit [km]
  const double ra_lag = 0.523599;     // Right ascension lag [rad]
  const int n_prm = 3;                // Harris-Priester parameter
                                      // 2(6) low(high) inclination

  // Variables

  int h_index = 0;                             // Height section variables
  double height;                               // Earth flattening
  double dec_Sun, ra_Sun, c_dec;               // Sun declination, right asc.
  double c_psi2;                               // Harris-Priester modification
  double density, h_min, h_max, d_min, d_max;  // Height, density parameters
  Vector r_Sun = Sun(mjd_TT);                  // Sun position
  Vector u(3);                                 // Apex of diurnal bulge

  // Satellite height
  height = Height(r_NP) / 1000.0;  //  [km]

  // Exit with zero density outside height model limits
  if (height >= upper_limit || height <= lower_limit) return 0.;

  // Sun right ascension, declination
  ra_Sun = atan2(r_Sun(1), r_Sun(0));
  dec_Sun = atan2(r_Sun(2), sqrt(r_Sun(0) * r_Sun(0) + r_Sun(1) * r_Sun(1)));

  // Unit vector u towards the apex of the diurnal bulge in inertial geocentric coordinates
  c_dec = cos(dec_Sun);
  u(0) = c_dec * cos(ra_Sun + ra_lag);
  u(1) = c_dec * sin(ra_Sun + ra_lag);
  u(2) = sin(dec_Sun);

  // Cosine of half angle between satellite position vector and apex of diurnal bulge
  c_psi2 = 0.5 + 0.5 * r_NP.dot(u) / r_NP.norm();

  // Height index search and exponential density interpolation
  for (int i = 0; i < N_Coeff_rho - 1; i++) {  // loop over N_Coeff_rho height regimes
    if (height >= Data_h[i] && height < Data_h[i + 1]) {
      h_index = i;  // h_index identifies height section
      break;
    }
  }

  h_min = (Data_h[h_index] - Data_h[h_index + 1]) / log(Data_rho_min[h_index + 1] / Data_rho_min[h_index]);
  h_max = (Data_h[h_index] - Data_h[h_index + 1]) / log(Data_rho_max[h_index + 1] / Data_rho_max[h_index]);

  d_min = Data_rho_min[h_index] * exp((Data_h[h_index] - height) / h_min);
  d_max = Data_rho_max[h_index] * exp((Data_h[h_index] - height) / h_max);

  // Density computation
  density = d_min + (d_max - d_min) * pow(c_psi2, n_prm);

  return density * 1.0e-12;  // [kg/m^3]
}

Vector AccelHarmonic(const Vector& r, const Matrix& E, int n_max, int m_max) {
  double r_sqr, rho, Fac;  // Auxiliary quantities
  double x0, y0, z0;       // Normalized coordinates
  double ax, ay, az;       // Acceleration vector
  double C, S;             // Gravitational coefficients
  Vector r_bf(3);          // Body-fixed position
  Vector a_bf(3);          // Body-fixed acceleration

  Matrix V(n_max + 2, n_max + 2);  // Harmonic functions
  Matrix W(n_max + 2, n_max + 2);  // work array (0..n_max+1,0..n_max+1)

  // Body-fixed position
  r_bf = E * r;
  // static int count = 0;
  // if (count == 0) {
  //   std::cout << "E = " << E << std::endl;
  //   std::cout << "r: " << r << std::endl;
  //   std::cout << "r_bf: " << r_bf << std::endl;
  //   count++;
  // }
  // Auxiliary quantities
  r_sqr = r_bf.dot(r_bf);  // Square of distance
  rho = R_Earth * R_Earth / r_sqr;

  // Normalized coordinates
  x0 = R_Earth * r_bf(0) / r_sqr;
  y0 = R_Earth * r_bf(1) / r_sqr;
  z0 = R_Earth * r_bf(2) / r_sqr;

  // Evaluate harmonic functions
  //   V_nm = (R_Earth/r)^(n+1) * P_nm(sin(phi)) * cos(m*lambda)
  // and
  //   W_nm = (R_Earth/r)^(n+1) * P_nm(sin(phi)) * sin(m*lambda)
  // up to degree and order n_max+1

  // Calculate zonal terms V(n,0); set W(n,0)=0.0
  V(0, 0) = R_Earth / sqrt(r_sqr);
  W(0, 0) = 0.0;

  V(1, 0) = z0 * V(0, 0);
  W(1, 0) = 0.0;

  for (int n = 2; n <= n_max + 1; n++) {
    V(n, 0) = ((2 * n - 1) * z0 * V(n - 1, 0) - (n - 1) * rho * V(n - 2, 0)) / n;
    W(n, 0) = 0.0;
  };

  // Calculate tesseral and sectorial terms
  for (int m = 1; m <= m_max + 1; m++) {
    // Calculate V(m,m) .. V(n_max+1,m)
    V(m, m) = (2 * m - 1) * (x0 * V(m - 1, m - 1) - y0 * W(m - 1, m - 1));
    W(m, m) = (2 * m - 1) * (x0 * W(m - 1, m - 1) + y0 * V(m - 1, m - 1));
    if (m <= n_max) {
      V(m + 1, m) = (2 * m + 1) * z0 * V(m, m);
      W(m + 1, m) = (2 * m + 1) * z0 * W(m, m);
    };
    for (int n = m + 2; n <= n_max + 1; n++) {
      V(n, m) = ((2 * n - 1) * z0 * V(n - 1, m) - (n + m - 1) * rho * V(n - 2, m)) / (n - m);
      W(n, m) = ((2 * n - 1) * z0 * W(n - 1, m) - (n + m - 1) * rho * W(n - 2, m)) / (n - m);
    };
  };

  // Calculate accelerations ax, ay, az (Cunningham recursion to compute the gradient of the potential)
  ax = ay = az = 0.0;
  for (int m = 0; m <= m_max; m++)
    for (int n = m; n <= n_max; n++)
      if (m == 0) {
        C = CS_JGM3(n, 0);  // = C_n,0
        ax -= C * V(n + 1, 1);
        ay -= C * W(n + 1, 1);
        az -= (n + 1) * C * V(n + 1, 0);
      } else {
        C = CS_JGM3(n, m);      // = C_n,m
        S = CS_JGM3(m - 1, n);  // = S_n,m
        Fac = 0.5 * (n - m + 1) * (n - m + 2);
        ax += +0.5 * (-C * V(n + 1, m + 1) - S * W(n + 1, m + 1)) + Fac * (+C * V(n + 1, m - 1) + S * W(n + 1, m - 1));
        ay += +0.5 * (-C * W(n + 1, m + 1) + S * V(n + 1, m + 1)) + Fac * (-C * W(n + 1, m - 1) + S * V(n + 1, m - 1));
        az += (n - m + 1) * (-C * V(n + 1, m) - S * W(n + 1, m));
      };

  // Body-fixed acceleration
  a_bf = Vector(ax, ay, az) * (GM_EARTH / (R_Earth * R_Earth));

  // Inertial acceleration
  return E.transpose() * a_bf;
}

int gravField(int n, double t, double x[], double f[], void* param) {
  if (n != 6) return 1;
  for (int i = 0; i < 3; i++) f[i] = x[i + 3];  // dx/dt = v_x, dy/dt = v_y, dz/dt = v_z
  double mjd_tt = t / 86400.0;
  args_gravField* prm = (args_gravField*)param;
  Vector r = Vector(x[0], x[1], x[2]);
  Vector v = Vector(x[3], x[4], x[5]);
  Vector v_dot = Vector(0, 0, 0);
  if (prm->pointEarth)
    v_dot += AccelPointEarth(r);
  else
    v_dot += AccelHarmonic(r, J20002ECEF(mjd_tt), prm->n_max, prm->m_max);
  Vector r_sun = Sun(mjd_tt);
  Vector v_dot_0 = v_dot;
  static int i = 0;
  // if (i < 100)
  //   std::cout << "v_dot_0 = " << v_dot.norm() << std::endl;
  if (prm->sun) v_dot += AccelPointMass(r, r_sun, GM_SUN);
  // if (i < 1000)
  //   std::cout << "v_dot_sun (" << prm->sun << ") = " << (v_dot - v_dot_0).norm() << std::endl;
  // v_dot_0 = v_dot;
  if (prm->moon) v_dot += AccelPointMass(r, Moon(mjd_tt), GM_MOON);
  // if (i < 1000)
  //   std::cout << "v_dot_moon (" << prm->moon << ") = " << (v_dot - v_dot_0).norm() << std::endl;
  v_dot_0 = v_dot;
  if (prm->otherPlanets) v_dot += AccelPointMass(r, Mars(mjd_tt), GM_MARS) + AccelPointMass(r, Venus(mjd_tt), GM_VENUS) +
                                  AccelPointMass(r, Jupiter(mjd_tt), GM_JUPITER);
  // if (i < 1000)
  //   std::cout << "v_dot_otherPlanets (" << prm->otherPlanets << ") = " << (v_dot - v_dot_0).norm() << std::endl;
  if (prm->solar_rad) v_dot += AccelSolarRad(r, r_sun, prm->Am) * Illumination(r, r_sun);
  // if (i < 100)
  //   std::cout << "v_dot_solar_rad (" << prm->solar_rad << ") = " << (v_dot - v_dot_0).norm() << std::endl;

  if (prm->atmo_drag) v_dot += AccelDrag(r, v, mjd_tt, NutationMatrix(mjd_tt) * PrecessionMatrix(mjd_tt), prm->Am);
  // if (i < 100)
  // std::cout << "v_dot_drag (" << prm->atmo_drag << ") = " << (v_dot - v_dot_0).norm() << std::endl;
  i++;
  // static int count = 0;
  // // std::cout << "count = " << count << std::endl;
  // if (count == 0) {
  //   std::cout << "v = " << Vector(x[3], x[4], x[5]) << std::endl;
  //   std::cout << "v_dot = " << v_dot << std::endl;
  //   count++;
  // }
  // if (prm->sun) v_dot += AccelPointMass(r,  GM_SUN);
  // if (prm->moon) v_dot += AccelPointMass(r, GM_MOON);
  for (int i = 0; i < 3; i++) f[i + 3] = v_dot(i);
  return 0;
}

int integrateOrbit(Satellite& s, double T, double h, double tol, int maxNumSteps, void* prm) {
  int n = 6;
  double x[6] = {s.r_ECI(0), s.r_ECI(1), s.r_ECI(2), s.v_ECI(0), s.v_ECI(1), s.v_ECI(2)};
  double hmin = h * 0.05, hmax = h * 20.0;

  double t = s.mjd_TT * 86400.0;
  // std::cout << "r_ECI = " << s.r_ECI << std::endl;
  // std::cout << "v_ECI = " << s.v_ECI << std::endl;
  // set precision
  // std::cout.precision(16);
  // std::cout << "t (mjd_tt) = " << t / 86400 << std::endl;
  if (flow(&t, x, &h, T, hmin, hmax, tol, maxNumSteps, n, gravField, prm)) return 1;

  s.mjd_TT = t / 86400.0;
  s.r_ECI = Vector(x[0], x[1], x[2]);
  s.v_ECI = Vector(x[3], x[4], x[5]);

  // s.set_orbital_elements(s.r_ECI, s.v_ECI);
  return 0;
}