#ifndef KEPLER_H
#define KEPLER_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  double e;  // eccentricity
  double M;  // mean anomaly
} args_kepler;

double kepler_equation(double E, void *args);
double diff_kepler_equation(double E, void *args);

#ifdef __cplusplus
}
#endif
#endif  // KEPLER_H