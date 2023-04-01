#include "TT_to_TDB.h"

#include "monte/SAT_Const.h"
double TT_to_TDB(double jd_TT_large, double jd_TT_small) {
  /*
******************************************
JD_TDB = JD_TT + P

// P is a periodic correction.
******************************************
*/
  double jd_TDB, P = 0;  // P is a periodic correction.
  double T;              // Julian centuries in TT
#define numCoeffs 5
  const double C[numCoeffs][3] =
      {
          // terms of the form A * sin(B * T + C)
          // A are in seconds
          // B, C are in degrees
          // A         B          C
          //
          {0.0016568, 35999.37, 357.5},  //   1
          {0.0000224, 32964.5, 246},     //   2
          {0.0000138, 71998.7, 355},     //   4
          {0.0000048, 3034.9, 25},       //   3
          {0.0000047, 34777.3, 230},     //   5
      };

  T = ((jd_TT_large - 2451545) + jd_TT_small) * 1. / 36525;
  for (int i = 0; i < numCoeffs; i++) {
    P += C[i][0] / Arcs * sin((C[i][1] * T + C[i][2]) / Deg);
  }
  // jd_TDB = jd_TT_large + jd_TT_small + P;

  // return jd_TDB;
  return P;
}
