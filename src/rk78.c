#include "../include/rk78.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

int rk78(double *t, double x[], double *h,
         double hmin, double hmax, double tol,
         int n, int (*field)(int n, double t, double x[], double f[], void *param),
         void *param) {
  // Local variables
  int return_value;  // to store the return value of field
  int indexA;        // to store the index of A
  double tt, dif, RK8norm, eps, temp;
  double K[13 * n];  // matrix vector containing all the stages
  double RK7[n];     // vector containing the RK7 step
  double RK8[n];     // vector containing the RK8 step

  do {  // do a step for both the RK7 and RK8 methods
    indexA = 0;
    // loop over stages (RK7 has 11 stages and RK8 has 13 stages)
    for (int i = 0; i < 13; i++) {
      memcpy(RK7, x, n * sizeof(double));  // we use RK7 as a temporary vector (in order to avoid using more memory)
      tt = (*t) + c[i] * (*h);             // time at stage i
      for (int j = 0; j < i; j++, indexA++) {
        temp = A[indexA] * (*h);     // temp = step size * A[i][j]
        for (int m = 0; m < n; m++)  // sum over all the components of the contribution of stage j
          RK7[m] += temp * K[j * n + m];
      }
      return_value = field(n, tt, RK7, K + i * n, param);
      if (return_value) return return_value;  // Error
    }
    RK8norm = dif = 0;
    memcpy(RK8, x, n * sizeof(double));
    memcpy(RK7, x, n * sizeof(double));

    // we loop over the coordinates and then over the stages to save time and compute d and dd at the same time.
    for (int m = 0; m < n; m++) {     // loop over the coordinates
      for (int i = 0; i < 11; i++) {  // loop over the first 11 stages that both RK7 and RK8 'share'
        temp = (*h) * K[i * n + m];
        RK7[m] += temp * brk7[i];
        RK8[m] += temp * brk8[i];
      }
      // add the contribution of the last two stages (on RK8)
      RK8[m] += (*h) * (brk8[11] * K[11 * n + m] + brk8[12] * K[12 * n + m]);
      dif += fabs(RK8[m] - RK7[m]);
      RK8norm += fabs(RK8[m]);  // L1 norm of RK8 used to relativize the tolerance
    }
    dif /= n;                          // average of the differences between RK7 and RK8
    eps = tol * (1. + .01 * RK8norm);  // updated tolerance that accounts for large values of the field (in that case the tolerance is increased).
    // if the error is small enough (less than our 'artifical' tolerance) or the step size is smaller than the minimum step size we break the loop
    // printf("dif = %g, eps = %g, *h = %g, tol = %g, RK8norm = %g\n", dif, eps, *h, tol, RK8norm);
    if (dif < eps || fabs(*h) <= fabs(hmin)) break;
    (*h) *= .9 * pow(eps / dif, .125);  // Fehlberg correction. Note that 0.125 = 1/(1 + 7). See Stoer 3rd edition 7.2.5.13. We take the conservative correction coefficient as 0.9.
    if (fabs(*h) < fabs(hmin)) *h = hmin;
  } while (1);
  (*t) += (*h);                        // we update the time
  memcpy(x, RK8, n * sizeof(double));  // we update the vector x with the RK8 step
  // If the error is super small, we don't want the next step to be huge. So we bound the error from below by eps / 256.
  if (dif < eps / 256) dif = eps / 256;

  (*h) *= .9 * pow(eps / dif, .125);  // Fehlberg correction. Note that 0.125 = 1/(1 + 7). See Stoer 3rd edition 7.2.5.13. We take the conservative correction coefficient as 0.9.
  // we impose the limits on the step size
  if (fabs(*h) < fabs(hmin))
    *h = hmin;
  else if (fabs(*h) > fabs(hmax))
    *h = hmax;
  return 0;
}

int flow(double *t, double x[], double *h, double T, double hmin, double hmax, double tol, int maxNumSteps, int n, int (*field)(int n, double t, double x[], double f[], void *param), void *param) {
  double t0 = *t;
  int count = 0;

  if (fabs(T) < tol) {
    *t = t0 + T;
    return 0;
  }
  printf("------------------------------\n--------------------------\n-----------------------\n");
  printf("T = %g, t0 = %g, *t = %g, h = %g, hmin = %g, hmax = %g \n", T, t0, *t, *h, hmin, hmax);
  while (fabs(*t - t0) < fabs(T) && count < maxNumSteps) {
    if (fabs(*t + *h - t0) > fabs(T)) {
      *h = t0 + T - *t;
      if (rk78(t, x, h, hmin, hmax, tol, n, field, param))
        return 1;
      else
        break;
    }
    printf("Hola\n");
    if (rk78(t, x, h, hmin, hmax, tol, n, field, param)) return 1;
    count++;
    printf("t = %g, t - t0 = %g, *t + *h - t0 = %g, T = %g, h = %g, hmin = %g, hmax = %g \n", *t, *t - t0, *t + *h - t0, T, *h, hmin, hmax);
  }
  printf("COOOOOUNT = %d\n", count);
  if (count == maxNumSteps) {
    return 1;
  } else {
    *t = t0 + T;
    return 0;
  }
}
