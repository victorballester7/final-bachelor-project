#include "../include/ReferenceSystems.h"

#include "../include/LinearAlgebra.h"
#include "../include/constants.h"

Matrix Perifocal2ECI(double i, double Omega, double omega) {
  Matrix R(3, -Omega), S(1, -i), T(3, -omega);
  return R * S * T;
}

Matrix J20002ECEF(double mjd_TT) {
  // double mjd_UT1 = mjd_TT;  // we assume that UT1 = TT
  static int count = 0;
  if (count == 0) {
    count++;
    // std::cout << "P = " << PrecessionMatrix(mjd_TT) << std::endl;
    // std::cout << "N = " << NutationMatrix(mjd_TT) << std::endl;
    // std::cout << "R = " << RotationMatrix(mjd_TT) << std::endl;
    // std::cout << "Polar = " << PolarMotionMatrix(mjd_TT) << std::endl;
  }
  return PolarMotionMatrix(mjd_TT) * RotationMatrix(mjd_TT) * NutationMatrix(mjd_TT) * PrecessionMatrix(mjd_TT);
}

Matrix ECI2TEME(double mjd_TT) {
  return EquinoxMatrix(mjd_TT) * NutationMatrix(mjd_TT) * PrecessionMatrix(mjd_TT);
}

// double TT2JDcenturies_TT(double t) {
//   // TT (Terrestrial Time) Julian Date
//   return Mjd_TT + 2400000.5;
// }

double UTC2UT1(double mjd_utc) {  // valid from 01/01/2023 onwards.
  int day_before = floor(mjd_utc);
  // we will do a linear interpolation between the two values of UT1-UTC
  double f1 = eop[day_before][5];      // first value of UT1-UTC
  double f2 = eop[day_before + 1][5];  // second value of UT1-UTC
  double ut1 = (mjd_utc - day_before) * f2 - (mjd_utc - day_before - 1) * f1;
  return ut1;
}

double TT2UTC(double mjd_tt) {
  double sec_tt_utc = 32.184 + 37.0;  // valid from 01/01/2017 onwards.
  double mjd_utc = mjd_tt - sec_tt_utc / 86400.0;
  return mjd_utc;
}

double UTC2TT(double mjd_utc) {
  double sec_tt_utc = 32.184 + 37.0;  // valid from 01/01/2017 onwards.
  double mjd_tt = mjd_utc + sec_tt_utc / 86400.0;
  return mjd_tt;
}

double MeanObliquity(double mjd_TT) {
  const double T = (mjd_TT - MJD_J2000) / 36525.0;

  return DEG_TO_RAD * (23.43929111 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T / DEG_TO_ARCSEC);
}

double EqEquinoxes(double mjd_TT) {
  double dpsi, deps;  // Nutation angles

  NutationAngles(mjd_TT, dpsi, deps);
  // static int count = 0;
  // if (count == 0) {
  //   count++;
  //   std::cout << "dpsi = " << dpsi / DEG_TO_RAD << std::endl;
  //   std::cout << "deps = " << deps / DEG_TO_RAD << std::endl;
  //   std::cout << "meanObli = " << MeanObliquity(mjd_TT) / DEG_TO_RAD << std::endl;
  //   std::cout << "MeanObliquity(mjd_TT) + deps = " << (MeanObliquity(mjd_TT) + deps) / DEG_TO_RAD << std::endl;
  // }

  // Equation of the equinoxes
  const double T = (mjd_TT - MJD_J2000) / 36525.0;
  const double T2 = T * T;
  const double T3 = T2 * T;
  const double T4 = T3 * T;
  const double rev = 360.0 * 3600.0;  // arcsec/revolution

  // From Vallado:
  //   Om  mean longitude of the ascending node
  double Om = fmod(450160.398036 - 6962890.5431 * T + 7.4722 * T2 + 0.007702 * T3 - 0.00005939 * T4, rev) * ARCSEC_TO_RAD;
  double moon_perturbation = (0.00264 * sin(Om) + 0.000063 * sin(2.0 * Om)) * ARCSEC_TO_RAD;
  double eq_equinox = dpsi * cos(MeanObliquity(mjd_TT));
  eq_equinox += moon_perturbation;

  // double eq_equinox = dpsi * cos(MeanObliquity(mjd_TT)+deps);
  // eq_equinox += 0;

  return eq_equinox;
};

double GMST(double mjd_UT1) {
  // Variables
  double T, gmst;
  // double T_0;

  // Mean Sidereal Time
  // mjd_UT1_0 = floor(mjd_UT1);
  // UT1 = 86400.0 * (mjd_UT1 - mjd_UT1_0);  // seconds in UT1 of that day [s]
  // T_0 = (mjd_UT1_0 - MJD_J2000) / 36525.0;
  T = (mjd_UT1 - MJD_J2000) / 36525.0;

  // gmst = 24110.54841 + 8640184.812866 * T_0 + 1.002737909350795 * UT1 + (0.093104 - 6.2e-6 * T) * T * T;  // [s] - Montenbruck (bad)
  // gmst = 24110.54841 + 8640184.812866 * T + (0.093104 - 6.2e-6 * T) * T * T + 1.002737909350795 * DEG_TO_ARCSEC * UT1;  // [s] - How I think it should be (https://gssc.esa.int/navipedia/index.php/CEP_to_ITRF) + Vallado

  // gmst = 100.4606184 + 36000.77005361 * T + 0.00038793 * T * T - 2.6e-8 * T * T * T + 1.002737909350795 * UT1;  // [deg] - Vallado

  gmst = 67310.54841 + (876600.0 * HOUR_TO_SEC + 8640184.812866) * T + 0.093104 * T * T - 6.2e-6 * T * T * T;  // [s] - Vallado

  // static int count = 0;
  // if (count == 0) {
  //   count++;
  //   // change precision output
  //   std::cout.precision(15);
  //   std::cout << "mjd_UT1 = " << mjd_UT1 << std::endl;
  //   std::cout << "mjd_UT1_0 = " << mjd_UT1_0 << std::endl;
  //   std::cout << "mjd_UT1 - mjd_UT1_0 = " << mjd_UT1 - mjd_UT1_0 << std::endl;
  //   std::cout << "UT1 = " << UT1 << std::endl;
  //   // std::cout << "T_0 = " << T_0 << std::endl;
  //   std::cout << "T = " << T << std::endl;
  //   std::cout << "gmst = " << gmst << std::endl;
  //   std::cout << "gmst_to_deg = " << fmod(gmst * SEC_TO_ARCSEC * ARCSEC_TO_RAD, 360.) << std::endl;
  // }
  // first we extract the fractional part of the gmst (which is the fraction of a day (or revolution)), and then we convert it to radians. The passing reference gmst is just to avoid a warning, but here we don't need it (it would be the integer part of gmst).
  // return modf(gmst / 86400.0, &gmst) * 2 * PI;  // [rad] - Correction based on Montenbruck
  // we sum 2*PI to make sure that the result is positive. Otherwise, the fmod function can also return a negative value.
  return fmod(fmod(gmst, 86400.) * SEC_TO_RAD + 2 * PI, 2 * PI);  // [rad] - Vallado

  // return fmod(gmst * SEC_TO_ARCSEC * ARCSEC_TO_RAD, 2 * PI);  // [rad] - My version (I checked that it is the same as the one above)
  // return fmod(gmst * DEG_TO_RAD, 2 * PI);  // [rad] - My version (I checked that it is the same as the one above)
}

double GAST(double mjd_TT) {
  // we assume UT1 = UTC
  double mjd_UT1 = TT2UTC(mjd_TT);
  mjd_UT1 = get_mjd(2004, 4, 6.3274067829);
  // static int count = 0;
  // if (count == 0) {
  //   count++;
  //   std::cout << "equinox = " << EqEquinoxes(mjd_TT) / DEG_TO_RAD << std::endl;
  //   std::cout << "gmstGAST = " << GMST(mjd_UT1) / DEG_TO_RAD << std::endl;
  //   std::cout << "mjd_TTGAST = " << mjd_TT << std::endl;
  // }
  return fmod(GMST(mjd_UT1) + EqEquinoxes(mjd_TT), 2 * PI);
}

void get_gregorian_date(int year, double day_frac, int& month, double& day) {
  int days_in_month[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  if (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0))
    days_in_month[1] = 29;  // Leap year

  day = day_frac;
  month = 0;
  while (day >= days_in_month[month]) {
    day -= days_in_month[month];
    month++;
  }
  month++;
}

double get_mjd(double yyyy, double mm, double dd) {
  // more info at
  // https://en.wikipedia.org/wiki/Julian_day#Converting_Julian_or_Gregorian_calendar_date_to_Julian_Day_Number
  // https://quasar.as.utexas.edu/BillInfo/JulianDatesG.html
  double a, b;
  if (mm <= 2) {
    yyyy--;
    mm += 12;
  }
  a = floor(yyyy / 100);
  b = yyyy + 0.01 * mm + 0.0001 * dd >= 1582.1015 ? 2 - a + floor(a / 4) : 0;  // Gregorian reform correction
  double year_part = floor(365.25 * (yyyy + 4716)) - 2400000.5;                // we do this in order not to lose precision
  return year_part + floor(30.6001 * (mm + 1)) + dd + b - 1524.5;
}

Matrix PrecessionMatrix(double mjd_TT) {
  // Constants
  const double T = (mjd_TT - MJD_J2000) / 36525.0;

  // std::cout << "T_prec = " << T << std::endl;

  // Variables
  double zeta, z, theta;

  // Precession angles
  zeta = (2306.2181 + (0.30188 + 0.017998 * T) * T) * T * ARCSEC_TO_RAD;
  theta = (2004.3109 - (0.42665 + 0.041833 * T) * T) * T * ARCSEC_TO_RAD;
  // z = zeta + (0.79280 + 0.000205 * T) * T * T * ARCSEC_TO_RAD;  // montenbruck
  z = (2306.2181 + (1.09468 + 0.018203 * T) * T) * T * ARCSEC_TO_RAD;  // vallado, I checked and it's the same as montenbruck

  // std::cout << "zeta = " << zeta / DEG_TO_RAD << std::endl;
  // std::cout << "z = " << z / DEG_TO_RAD << std::endl;
  // std::cout << "theta = " << theta / DEG_TO_RAD << std::endl;

  // static int count = 0;

  // if (count == 0) {
  // count++;
  // std::cout << "MJD_TT: " << mjd_TT << std::endl;
  // std::cout << "T: " << T << std::endl;
  // std::cout << "zeta: " << zeta << std::endl;
  // std::cout << "z: " << z << std::endl;
  // std::cout << "theta: " << theta << std::endl;
  // }

  // Precession matrix
  Matrix R(3, -z), S(2, theta), V(3, -zeta);
  return R * S * V;
}

void NutationAngles(double mjd_TT, double& dpsi, double& deps) {
  // Constants
  const double T = (mjd_TT - MJD_J2000) / 36525.0;
  const double T2 = T * T;
  const double T3 = T2 * T;
  const double T4 = T3 * T;
  const double rev = 360.0 * 3600.0;  // arcsec/revolution

  int nut_coeffs = 106;  // number of coefficients in nutation model (must be <= 106)

  // Variables
  double l, lp, F, D, Om;
  double arg;

  // Mean arguments of luni-solar motion
  //
  //   l   mean anomaly of the Moon
  //   l'  mean anomaly of the Sun
  //   F   mean argument of latitude
  //   D   mean longitude elongation of the Moon from the Sun
  //   Om  mean longitude of the ascending node

  // From Montebruck
  // l = fmod(485866.733 + (1325.0 * rev + 715922.633) * T + 31.310 * T2 + 0.064 * T3, rev);
  // lp = fmod(1287099.804 + (99.0 * rev + 1292581.224) * T - 0.577 * T2 - 0.012 * T3, rev);
  // F = fmod(335778.877 + (1342.0 * rev + 295263.137) * T - 13.257 * T2 + 0.011 * T3, rev);
  // D = fmod(1072261.307 + (1236.0 * rev + 1105601.328) * T - 6.891 * T2 + 0.019 * T3, rev);
  // Om = fmod(450160.280 - (5.0 * rev + 482890.539) * T + 7.455 * T2 + 0.008 * T3, rev);

  // From Vallado:
  l = fmod(485868.249036 + 1717915923.2178 * T + 31.8792 * T2 + 0.051635 * T3 - 0.00024470 * T4, rev);
  lp = fmod(1287104.79305 + 129596581.0481 * T - 0.5532 * T2 - 0.000136 * T3 - 0.00001149 * T4, rev);
  F = fmod(335779.526232 + 1739527262.8478 * T - 12.7512 * T2 - 0.001037 * T3 + 0.00000417 * T4, rev);
  D = fmod(1072260.70369 + 1602961601.2090 * T - 6.3706 * T2 + 0.006593 * T3 - 0.00003169 * T4, rev);
  Om = fmod(450160.398036 - 6962890.5431 * T + 7.4722 * T2 + 0.007702 * T3 - 0.00005939 * T4, rev);

  // Nutation in longitude and obliquity [rad]

  deps = dpsi = 0.0;
  for (int i = 0; i < nut_coeffs; i++) {
    arg = (C[i][0] * l + C[i][1] * lp + C[i][2] * F + C[i][3] * D + C[i][4] * Om) * ARCSEC_TO_RAD;
    dpsi += (C[i][5] + C[i][6] * T) * sin(arg);
    deps += (C[i][7] + C[i][8] * T) * cos(arg);
  };
  dpsi *= 1.0E-5;  // Now they are in arcseconds
  deps *= 1.0E-5;  // Now they are in arcseconds

  // We add the corrections of the IERS Bulletin B
  double mjd_UTC = TT2UTC(mjd_TT);
  int index = 0;
  for (int i = 0; i < eop_days; i++) {
    if (mjd_UTC >= eop[i][3])
      continue;
    else
      index = i - 1;
  }
  // we do a linear interpolation to get the 'exact' xp and yp at the given time

  double day = floor(mjd_UTC);
  double ddpsi = eop[index + 1][6] * (mjd_UTC - day) - eop[index][6] * (mjd_UTC - day - 1);
  double ddeps = eop[index + 1][7] * (mjd_UTC - day) - eop[index][7] * (mjd_UTC - day - 1);

  ddpsi *= 1.0E-3;  // Now they are in arcseconds
  ddeps *= 1.0E-3;  // Now they are in arcseconds

  dpsi += ddpsi;
  deps += ddeps;

  dpsi *= ARCSEC_TO_RAD;
  deps *= ARCSEC_TO_RAD;
}

Matrix NutationMatrix(double mjd_TT) {
  double dpsi, deps, eps;

  // Mean obliquity of the ecliptic
  eps = MeanObliquity(mjd_TT);

  // Nutation in longitude and obliquity
  NutationAngles(mjd_TT, dpsi, deps);

  // Transformation from mean to true equator and equinox
  Matrix R(1, -eps - deps), S(3, -dpsi), T(1, eps);
  return R * S * T;
}

Matrix EquinoxMatrix(double mjd_TT) {
  return Matrix(3, EqEquinoxes(mjd_TT));
}

Matrix RotationMatrix(double mjd_TT) {
  // we assume UT1 = UTC
  // double incr = (32.184 + 37.0) / 86400.0;  // TT-UT1 [days]
  // double mjd_UT1 = mjd_TT - incr;
  // // return Matrix(3, GMST(mjd_UT1)) * Matrix(3, EqEquinoxes(mjd_TT));
  return Matrix(3, GAST(mjd_TT));
}

// Matrix RotationMatrix(double mjd_UT1) {
//   return Matrix(3, GAST(mjd_UT1));
// }

void PolarMotion(double mjd_TT, double& xp, double& yp) {
  double mjd_UTC = TT2UTC(mjd_TT);
  int index = 0;
  for (int i = 0; i < eop_days; i++) {
    if (mjd_UTC >= eop[i][3])
      continue;
    else
      index = i - 1;
  }
  // we do a linear interpolation to get the 'exact' xp and yp at the given time

  double day = floor(mjd_UTC);
  xp = eop[index + 1][4] * (mjd_UTC - day) - eop[index][4] * (mjd_UTC - day - 1);
  yp = eop[index + 1][5] * (mjd_UTC - day) - eop[index][5] * (mjd_UTC - day - 1);

  xp *= ARCSEC_TO_RAD;
  yp *= ARCSEC_TO_RAD;
}

Matrix PolarMotionMatrix(double mjd_TT) {
  double xp, yp;

  // Polar motion
  PolarMotion(mjd_TT, xp, yp);

  // Transformation from ITRS to GCRS
  Matrix R(2, -xp), S(1, -yp);
  return R * S;
}

double Height(Vector r_PN) {
  const double eps = 10e-7;  // Convergence criterion

  const double X = r_PN(0);  // Cartesian coordinates
  const double Y = r_PN(1);
  const double Z = r_PN(2);
  const double rho2 = X * X + Y * Y;  // Square of distance from z-axis

  // Iteration for obtaining dZ

  double dZ, dZ_new, sinPhi;
  double ZdZ, Nh, N;

  dZ = e2 * Z;
  while (true) {
    ZdZ = Z + dZ;
    Nh = sqrt(rho2 + ZdZ * ZdZ);
    sinPhi = ZdZ / Nh;  // Sine of geodetic latitude
    N = R_Earth / sqrt(1.0 - e2 * sinPhi * sinPhi);
    dZ_new = N * e2 * sinPhi;
    if (fabs(dZ - dZ_new) < eps) break;
    dZ = dZ_new;
  }

  // Longitude, latitude, altitude
  // lon = atan2(Y, X);
  // lat = atan2(ZdZ, sqrt(rho2));
  // h = Nh - N;
  return Nh - N;
}

Vector Sun(double mjd_tt) {
  // we asume mjd_ut1 = mjd_tt

  const double eps = MeanObliquity(mjd_tt);         // Mean obliquity of the ecliptic
  const double T = (mjd_tt - MJD_J2000) / 36525.0;  // Julian cent. since J2000

  double lM, L, M, r;
  Vector r_Sun(3);

  // Mean anomaly, ecliptic longitude and radius (from Vallado, p. 304)
  lM = (280.460 + 36000.770 * T) * DEG_TO_RAD;                                 // [rad]
  M = (357.5291092 + 35999.05034 * T) * DEG_TO_RAD;                            // [rad]
  r = (1.000140612 - 0.016708617 * cos(M) - 0.000139589 * cos(2.0 * M)) * AU;  // [m]

  lM = fmod(lM, 2 * PI);
  M = fmod(M, 2 * PI);

  // Ecliptic longitude
  L = lM + (1.914666471 * sin(M) + 0.019994643 * sin(2.0 * M)) * DEG_TO_RAD;  // [rad]
  L = fmod(L, 2 * PI);

  static int count = 0;
  if (count == 0) {
    std::cout << "T: " << T << std::endl;
    std::cout << "eps: " << eps / DEG_TO_RAD << std::endl;
    std::cout << "lM: " << fmod(lM, 2 * PI) / DEG_TO_RAD << std::endl;
    std::cout << "M: " << fmod(M, 2 * PI) / DEG_TO_RAD << std::endl;
    std::cout << "r: " << r / AU << std::endl;
    std::cout << "L: " << fmod(L, 2 * PI) / DEG_TO_RAD << std::endl;
    std::cout << "hola" << std::endl;
    std::cout << "position: " << Matrix(1, -eps) * Vector(r * cos(L), r * sin(L), 0.0) << std::endl;
    count++;
  }
  // Equatorial position vector in Mean-of-Date frame
  r_Sun = Matrix(1, -eps) * Vector(r * cos(L), r * sin(L), 0.0);

  // Equatorial position vector in J2000 frame
  r_Sun = PrecessionMatrix(mjd_tt).transpose() * r_Sun;

  return r_Sun;
}

Vector Moon(double mjd_TT) {
  // we asume mjd_ut1 = mjd_tt

  const double eps = MeanObliquity(mjd_TT);         // Mean obliquity of the ecliptic
  const double T = (mjd_TT - MJD_J2000) / 36525.0;  // Julian cent. since J2000

  double lM, M_moon, M_sun, D_sun, um_moon, Lat, Long, P, R;
  Vector r_Moon(3);

  // Mean elements of lunar orbit (from Vallado, p. 210)
  lM = (218.32 + 481267.883 * T) * DEG_TO_RAD;                                                                            // [rad]
  M_moon = (485868.249036 + (1717915923.2178 + (31.8792 + (0.051635 - 0.00024470 * T) * T) * T) * T) * ARCSEC_TO_RAD;     // [rad]
  M_sun = (1287104.79305 + (129596581.0481 + (-0.5532 + (0.000136 - 0.00001149 * T) * T) * T) * T) * ARCSEC_TO_RAD;       // [rad]
  D_sun = (1072260.70369 + (1602961601.2090 + (-6.3706 + (0.006593 - 0.00003169 * T) * T) * T) * T) * ARCSEC_TO_RAD;      // [rad]
  um_moon = (335779.526232 + (1739527262.8478 + (-12.7512 + (-0.001037 + 0.00000417 * T) * T) * T) * T) * ARCSEC_TO_RAD;  // [rad]

  // Ecliptic longitude
  Long = lM + (6.288750 * sin(M_moon) + 1.274018 * sin(2 * D_sun - M_moon) + 0.658309 * sin(2 * D_sun) + 0.213616 * sin(2 * M_moon) - 0.185596 * sin(M_sun) - 0.114336 * sin(2 * um_moon)) * DEG_TO_RAD;  // [rad]

  // Ecliptic latitude
  Lat = (5.128189 * sin(D_sun) + 0.280606 * sin(M_moon + D_sun) + 0.277693 * sin(M_moon - D_sun) + 0.173238 * sin(2 * D_sun - M_moon)) * DEG_TO_RAD;  // [rad]

  // Parallax
  P = (0.9508 + 0.0518 * cos(M_moon) + 0.0095 * cos(2 * D_sun - M_moon) + 0.0078 * cos(2 * D_sun) + 0.0028 * cos(2 * M_moon)) * DEG_TO_RAD;  // [rad]

  // Distance [m]

  R = R_Earth / sin(P);

  // Equatorial coordinates Mean-of-date
  r_Moon = Matrix(1, -eps) * Vector(R * cos(Long) * cos(Lat), R * sin(Long) * cos(Lat), R * sin(Lat));

  // Equatorial coordinates J2000
  r_Moon = PrecessionMatrix(mjd_TT).transpose() * r_Moon;

  return r_Moon;
}
