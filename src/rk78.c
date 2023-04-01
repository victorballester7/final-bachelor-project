#include "../include/rk78.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

int rk78(double *t, double x[], double *h,
         double hmin, double hmax, double tol,
         int n, int (*field)(int n, double t, double x[], double f[], void *param),
         void *param) {
  // Local variables
  int ib, j, k, l, iret;
  double tt, bet, d, dd, e3, wksp[15 * n], *r = wksp, *b = wksp + 13 * n,
                                           *f = wksp + 14 * n;
  /* Bucle en el pas */
  do {
    /* Fem RK7 -> b[] i RK8 -> f[]
     * Promig per coordenades de RK7-RK8 -> d
     * Norma sub-1 de RK7-RK8 -> dd */
    ib = 0;
    for (j = 0; j < 13; j++) {
      memcpy(b, x, n * sizeof(double));
      tt = (*t) + alfa[j] * (*h);
      for (k = 0; k < j; k++, ib++) {
        bet = beta[ib] * (*h);
        for (l = 0; l < n; l++)
          b[l] += bet * r[n * k + l];
      }
      iret = field(n, tt, b /*x*/, r + n * j /*f*/, param);
      if (iret) return iret;
    }
    dd = d = 0;
    for (l = 0; l < n; l++) {
      b[l] = f[l] = x[l];
      for (k = 0; k < 11; k++) {
        bet = (*h) * r[k * n + l];
        b[l] += bet * c[k];
        f[l] += bet * cp[k];
      }
      f[l] += (*h) * (cp[11] * r[11 * n + l] + cp[12] * r[12 * n + l]);
      d += fabs(f[l] - b[l]);
      dd += fabs(f[l]);
    }
    d /= n;                     /* No es la norma sub-1, sino el promig dels errors	*/
                                /*    comesos a cada coordenada. dd si que es la norma	*/
                                /*    sub-1, pero no s'usa per mesurar errors sino per	*/
                                /*    relativitzar la tolerancia.			*/
    e3 = tol * (1. + .01 * dd); /* Es una tolerancia absoluta per valors    */
                                /* petits de les coordenades, relativa      */
                                /* "retardada dos digits" per valors grans. */
                                /* Pleguem si O.K. o pas minim. */
    if (d < e3 || fabs(*h) <= hmin) break;
    /* Corregim pas per a reintegrar. */
    (*h) *= .9 * pow(e3 / d, .125);
    if (fabs((*h)) < hmin) (*h) = SGN(*h) * hmin;
    /* Torno a fer RK78 */
  } while (1);
  /* Guardem temps final. */
  (*t) += (*h);
  /* Guardem punt final. */
  memcpy(x, f, n * sizeof(double));
  /* Fem correccio de pas */
  if (d < e3) {
    double tmp = e3 / 256;
    if (d < tmp) d = tmp;
  }
  /* Si l'error ha estat molt petit	  */
  /*    no volem que el nou pas es dispari. */
  (*h) *= .9 * pow(e3 / d, .125); /* Correccio Fehlberg (Stoer 2a ed (7.2.5.16), */
                                  /*   noti's que NO es la (7.2.5.17)).          */
  if (fabs(*h) < hmin) {          /* Fem que estigui dins */
    (*h) = SGN(*h) * hmin;
    // fprintf(stderr, "rk78():: t %G : ajusto a pasmin %G !!\n", *t, hmin);
  } else if (fabs(*h) > hmax) { /* els limits permesos. */
    (*h) = SGN(*h) * hmax;
    // fprintf(stderr, "rk78():: t %G : ajusto a pasmax %G\n", *t, hmax);
  }
  return 0;
}

int flow(double *t, double x[], double *h, double T, double hmin, double hmax, double tol, int maxNumSteps, int n, int (*field)(int n, double t, double x[], double f[], void *param), void *param) {
  double t0 = *t;
  int count = 0;
  while (*t < t0 + T && count < maxNumSteps) {
    if (*t + *h > t0 + T)
      *h = t0 + T - *t;
    rk78(t, x, h, hmin, hmax, tol, n, field, param);
    count++;
  }
  if (*t < t0 + T - tol) {  // this means count = maxNumSteps
    return -1;
  } else {
    *t = t0 + T;
    return 0;
  }
}
