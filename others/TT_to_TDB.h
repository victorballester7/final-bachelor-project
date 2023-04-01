#include <math.h>
#include <stdio.h>

/*
Converts Julian date in TT (Terrestial Time) to Julian date in TDB (Barycentric Dynamical Time).

Arguments:
  - jd_TT_large: big digits of the Julian date in TT
  - jd_TT_small: small digits of the Julian date in TT

  For example if we want to compute the date 2460010.25 we could do:
   jd_TT_large + jd_TT_small = 2459984.5 + 25.75 = 2460010.25
   (Note that 2459984.5 is exact in binary notation and so we do not lose precision we doing substractions of large numbers)

Return: Julian date in Barycentric Dynamical Time
*/
double TT_to_TDB(double jd_TT_large, double jd_TT_small);
