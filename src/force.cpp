#include "../include/force.h"

#include <cmath>

#include "../include/LinearAlgebra.h"
#include "../include/ReferenceSystems.h"
#include "../include/Satellite.h"
#include "../include/constants.h"
#include "../include/rk78.h"

Vector AccelPointMass(const Vector& r, const Vector& s, double GM) {
  // relative distance between big body and small body
  Vector r_rel = s - r;
  double r3 = pow(r_rel.norm(), 3.);
  double s3 = pow(s.norm(), 3.);
  Vector sat_pointMass = r_rel * (GM / r3);
  Vector earth_pointMass = s * (GM / s3);
  return sat_pointMass - earth_pointMass;
}

Vector AccelSolarRad(const Vector& r, const Vector& r_Sun, double Area, double mass) {
  Vector d(3);

  // Relative position vector of spacecraft w.r.t. Sun

  d = r - r_Sun;

  double d3 = pow(d.norm(), 3);  // Distance cubed

  // Acceleration

  return d * (C_R * (Area / mass) * P0 * AU * AU / d3);
}

double Illumination(const Vector& r, const Vector& r_Sun) {
  Vector e_Sun = r_Sun / r_Sun.norm();  // Sun direction unit vector
  double s = r.dot(e_Sun);              // Projection of s/c position

  return ((s > 0 || (r - e_Sun * s).norm() > R_JGM3) ? 1.0 : 0.0);
}

// Vector AccelDrag(const Vector& r, const Vector& v, double mjd_TT, const Matrix& T, double Area, double mass) {
// }

// double Density_HP(double Mjd_TT, const Vector& r_tod) {
//   // Constants

//   const double upper_limit = 1000.0;  // Upper height limit [km]
//   const double lower_limit = 100.0;   // Lower height limit [km]
//   const double ra_lag = 0.523599;     // Right ascension lag [rad]
//   const int n_prm = 3;                // Harris-Priester parameter
//                                       // 2(6) low(high) inclination

//   // Harris-Priester atmospheric density model parameters
//   // Height [km], minimum density, maximum density [gm/km^3]

//   const int N_Coef = 50;
//   const double Data_h[N_Coef] = {
//       100.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0,
//       210.0, 220.0, 230.0, 240.0, 250.0, 260.0, 270.0, 280.0, 290.0, 300.0,
//       320.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0,
//       520.0, 540.0, 560.0, 580.0, 600.0, 620.0, 640.0, 660.0, 680.0, 700.0,
//       720.0, 740.0, 760.0, 780.0, 800.0, 840.0, 880.0, 920.0, 960.0, 1000.0};
//   const double Data_c_min[N_Coef] = {
//       4.974e+05, 2.490e+04, 8.377e+03, 3.899e+03, 2.122e+03, 1.263e+03,
//       8.008e+02, 5.283e+02, 3.617e+02, 2.557e+02, 1.839e+02, 1.341e+02,
//       9.949e+01, 7.488e+01, 5.709e+01, 4.403e+01, 3.430e+01, 2.697e+01,
//       2.139e+01, 1.708e+01, 1.099e+01, 7.214e+00, 4.824e+00, 3.274e+00,
//       2.249e+00, 1.558e+00, 1.091e+00, 7.701e-01, 5.474e-01, 3.916e-01,
//       2.819e-01, 2.042e-01, 1.488e-01, 1.092e-01, 8.070e-02, 6.012e-02,
//       4.519e-02, 3.430e-02, 2.632e-02, 2.043e-02, 1.607e-02, 1.281e-02,
//       1.036e-02, 8.496e-03, 7.069e-03, 4.680e-03, 3.200e-03, 2.210e-03,
//       1.560e-03, 1.150e-03};
//   const double Data_c_max[N_Coef] = {
//       4.974e+05, 2.490e+04, 8.710e+03, 4.059e+03, 2.215e+03, 1.344e+03,
//       8.758e+02, 6.010e+02, 4.297e+02, 3.162e+02, 2.396e+02, 1.853e+02,
//       1.455e+02, 1.157e+02, 9.308e+01, 7.555e+01, 6.182e+01, 5.095e+01,
//       4.226e+01, 3.526e+01, 2.511e+01, 1.819e+01, 1.337e+01, 9.955e+00,
//       7.492e+00, 5.684e+00, 4.355e+00, 3.362e+00, 2.612e+00, 2.042e+00,
//       1.605e+00, 1.267e+00, 1.005e+00, 7.997e-01, 6.390e-01, 5.123e-01,
//       4.121e-01, 3.325e-01, 2.691e-01, 2.185e-01, 1.779e-01, 1.452e-01,
//       1.190e-01, 9.776e-02, 8.059e-02, 5.741e-02, 4.210e-02, 3.130e-02,
//       2.360e-02, 1.810e-02};

//   const Vector h(&Data_h[0], N_Coef);
//   const Vector c_min(&Data_c_min[0], N_Coef);
//   const Vector c_max(&Data_c_max[0], N_Coef);

//   // Variables

//   int i, ih;                                   // Height section variables
//   double height;                               // Earth flattening
//   double dec_Sun, ra_Sun, c_dec;               // Sun declination, right asc.
//   double c_psi2;                               // Harris-Priester modification
//   double density, h_min, h_max, d_min, d_max;  // Height, density parameters
//   Vector r_Sun(3);                             // Sun position
//   Vector u(3);                                 // Apex of diurnal bulge

//   // Satellite height

//   height = Geodetic(r_tod).h / 1000.0;  //  [km]

//   // Exit with zero density outside height model limits

//   if (height >= upper_limit || height <= lower_limit) {
//     return 0.0;
//   }

//   // Sun right ascension, declination

//   r_Sun = Sun(Mjd_TT);
//   ra_Sun = atan2(r_Sun(1), r_Sun(0));
//   dec_Sun = atan2(r_Sun(2), sqrt(pow(r_Sun(0), 2) + pow(r_Sun(1), 2)));

//   // Unit vector u towards the apex of the diurnal bulge
//   // in inertial geocentric coordinates

//   c_dec = cos(dec_Sun);
//   u(0) = c_dec * cos(ra_Sun + ra_lag);
//   u(1) = c_dec * sin(ra_Sun + ra_lag);
//   u(2) = sin(dec_Sun);

//   // Cosine of half angle between satellite position vector and
//   // apex of diurnal bulge

//   c_psi2 = 0.5 + 0.5 * Dot(r_tod, u) / Norm(r_tod);

//   // Height index search and exponential density interpolation

//   ih = 0;                           // section index reset
//   for (i = 0; i < N_Coef - 1; i++)  // loop over N_Coef height regimes
//   {
//     if (height >= h(i) && height < h(i + 1)) {
//       ih = i;  // ih identifies height section
//       break;
//     }
//   }

//   h_min = (h(ih) - h(ih + 1)) / log(c_min(ih + 1) / c_min(ih));
//   h_max = (h(ih) - h(ih + 1)) / log(c_max(ih + 1) / c_max(ih));

//   d_min = c_min(ih) * exp((h(ih) - height) / h_min);
//   d_max = c_max(ih) * exp((h(ih) - height) / h_max);

//   // Density computation

//   density = d_min + (d_max - d_min) * pow(c_psi2, n_prm);

//   return density * 1.0e-12;  // [kg/m^3]
// }

Vector AccelHarmonic(const Vector& r, const Matrix& E, double GM, double R_ref, const Matrix& CS, int n_max, int m_max) {
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
  rho = R_ref * R_ref / r_sqr;

  // Normalized coordinates
  x0 = R_ref * r_bf(0) / r_sqr;
  y0 = R_ref * r_bf(1) / r_sqr;
  z0 = R_ref * r_bf(2) / r_sqr;

  // Evaluate harmonic functions
  //   V_nm = (R_ref/r)^(n+1) * P_nm(sin(phi)) * cos(m*lambda)
  // and
  //   W_nm = (R_ref/r)^(n+1) * P_nm(sin(phi)) * sin(m*lambda)
  // up to degree and order n_max+1

  // Calculate zonal terms V(n,0); set W(n,0)=0.0
  V(0, 0) = R_ref / sqrt(r_sqr);
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
        C = CS(n, 0);  // = C_n,0
        ax -= C * V(n + 1, 1);
        ay -= C * W(n + 1, 1);
        az -= (n + 1) * C * V(n + 1, 0);
      } else {
        C = CS(n, m);      // = C_n,m
        S = CS(m - 1, n);  // = S_n,m
        Fac = 0.5 * (n - m + 1) * (n - m + 2);
        ax += +0.5 * (-C * V(n + 1, m + 1) - S * W(n + 1, m + 1)) + Fac * (+C * V(n + 1, m - 1) + S * W(n + 1, m - 1));
        ay += +0.5 * (-C * W(n + 1, m + 1) + S * V(n + 1, m + 1)) + Fac * (-C * W(n + 1, m - 1) + S * V(n + 1, m - 1));
        az += (n - m + 1) * (-C * V(n + 1, m) - S * W(n + 1, m));
      };

  // Body-fixed acceleration
  a_bf = Vector(ax, ay, az) * (GM / (R_ref * R_ref));

  // Inertial acceleration
  return E.transpose() * a_bf;
}

int gravField(int n, double t, double x[], double f[], void* param) {
  if (n != 6) return 1;
  for (int i = 0; i < 3; i++) f[i] = x[i + 3];  // dx/dt = v_x, dy/dt = v_y, dz/dt = v_z
  double mjd_tt = t / 86400.0;
  args_gravField* prm = (args_gravField*)param;
  Vector r = Vector(x[0], x[1], x[2]);
  Vector v_dot = Vector(0, 0, 0);
  if (prm->pointEarth)
    v_dot += AccelPointMass(r, Vector(0, 0, 0), GM_EARTH);
  else
    v_dot += AccelHarmonic(r, J20002ECEF(mjd_tt), GM_EARTH, R_JGM3, CS_JGM3, prm->n_max, prm->m_max);
  Vector r_sun = Sun(mjd_tt);
  Vector r_moon = Moon(mjd_tt);
  Vector v_dot_0 = v_dot;
  // std::cout << "v_dot_0 = " << v_dot.norm() << std::endl;
  if (prm->sun) v_dot += AccelPointMass(r, r_sun, GM_SUN);
  // std::cout << "v_dot_sun = " << (v_dot - v_dot_0).norm() << std::endl;
  v_dot_0 = v_dot;
  if (prm->moon) v_dot += AccelPointMass(r, r_moon, GM_MOON);
  // std::cout << "v_dot_moon = " << (v_dot - v_dot_0).norm() << std::endl;
  v_dot_0 = v_dot;
  if (prm->solar_rad) v_dot += AccelSolarRad(r, r_sun, prm->A, prm->m) * Illumination(r, r_sun);
  // std::cout << "v_dot_solar_rad = " << (v_dot - v_dot_0).norm() << std::endl;

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
  double hmin = h * 0.5, hmax = h * 2.0;

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