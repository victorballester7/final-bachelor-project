#include <math.h>
#include <stdio.h>

int newton(double* x, double tol, int maxIter, double (*f)(double, void* param), double (*df)(double, void* param), void* param) {
  int iter = 0;
  double fx = f(*x, param);
  double dfx = df(*x, param);
  while (fabs(fx) > tol && iter < maxIter) {
    *x = *x - fx / dfx;
    fx = f(*x, param);
    dfx = df(*x, param);
    iter++;
  }
  if (iter == maxIter) {
    // printf("Maximum number of iterations exceeded.\n");
    return 1;
  }

  return 0;
}
