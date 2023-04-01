#include "../include/force.h"

#include <cmath>

#include "../include/LinearAlgebra.h"
#include "../include/ReferenceSystems.h"
#include "../include/Satellite.h"
#include "../include/constants.h"
#include "../include/rk78.h"

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
  for (int i = 0; i < 3; i++) f[i] = x[i + 3];
  double* prm = (double*)param;
  double n_max = prm[0];
  double m_max = prm[1];

  Vector v_dot = AccelHarmonic(Vector(x[0], x[1], x[2]), GCRF2ITRF(t), GM_EARTH, R_JGM3, CS_JGM3, n_max, m_max);
  for (int i = 0; i < 3; i++) f[i + 3] = v_dot(i);
  return 0;
}

int integrateOrbit(Satellite& s0, Satellite& sT, double T, double hmin, double hmax, double tol, int maxNumSteps, int n_max, int m_max) {
  args_gravField prm = {n_max, m_max};
  int n = 6;
  double x[6] = {s0.r_GCRF(0), s0.r_GCRF(1), s0.r_GCRF(2), s0.v_GCRF(0), s0.v_GCRF(1), s0.v_GCRF(2)};
  double h = 0.01 * SGN(T);

  printf("\n\nInitial position and velocity:\n");
  for (int i = 0; i < n; i++)
    printf("%lf ", x[i]);

  if (flow(&s0.mjd_TT, x, &h, T, hmin, hmax, tol, maxNumSteps, n, gravField, &prm)) return 1;

  printf("Final position and velocity:\n");
  for (int i = 0; i < n; i++)
    printf("%lf ", x[i]);

  sT.r_GCRF = Vector(x[0], x[1], x[2]);
  sT.v_GCRF = Vector(x[3], x[4], x[5]);
  sT.set_orbital_elements(sT.r_GCRF, sT.v_GCRF);
  return 0;
}
