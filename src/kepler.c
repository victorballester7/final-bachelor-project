#include "../include/kepler.h"

#include <math.h>
#include <stdio.h>

double kepler_equation(double E, void *args) {
  double *param = (double *)args;
  double e = param[0];
  double M = param[1];
  return E - e * sin(E) - M;
}

double diff_kepler_equation(double E, void *args) {
  double e = *((double *)args);
  return 1 - e * cos(E);
}
