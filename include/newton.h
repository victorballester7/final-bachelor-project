#ifndef NEWTON_H
#define NEWTON_H

#ifdef __cplusplus
extern "C" {
#endif
#include <math.h>
#include <stdarg.h>

// ----------------------------------------------
// Newton's method
// ----------------------------------------------
// x: initial guess and is where the solution is stored
// tol: tolerance
// maxIter: maximum number of iterations
// f: function to be solved
// df: derivative of the function
// param: parameters for the function (usually a pointer to a struct)
// returns: 0 if converged, 1 if not converged
// ----------------------------------------------
int newton(double* x, double tol, int maxIter, double (*f)(double, void* param), double (*df)(double, void* param), void* param);

#ifdef __cplusplus
}
#endif

#endif  // NEWTON_H
