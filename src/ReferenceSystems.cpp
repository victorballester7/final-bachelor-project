#include "../include/ReferenceSystems.h"

#include "../include/LinearAlgebra.h"
#include "../include/constants.h"

Matrix Perifocal2ECI(double i, double Omega, double omega) {
  Matrix R(3, -Omega), S(1, -i), T(3, -omega);
  return R * S * T;
}

Matrix ECI2ECEF(double mjd_TT) {
  // double mjd_UT1 = mjd_TT;  // we assume that UT1 = TT
  static int count = 0;
  if (count == 0) {
    count++;
    std::cout << "P = " << PrecessionMatrix(mjd_TT) << std::endl;
    std::cout << "N = " << NutationMatrix(mjd_TT) << std::endl;
    std::cout << "R = " << RotationMatrix(mjd_TT) << std::endl;
    std::cout << "Polar = " << PolarMotionMatrix(mjd_TT) << std::endl;
  }
  return PolarMotionMatrix(mjd_TT) * RotationMatrix(mjd_TT) * NutationMatrix(mjd_TT) * PrecessionMatrix(mjd_TT);
}

// double TT2JDcenturies_TT(double t) {
//   // TT (Terrestrial Time) Julian Date
//   return Mjd_TT + 2400000.5;
// }

double UTC2UT1(double utc) {
  int day_before = floor(utc);
  // we will do a linear interpolation between the two values of UT1-UTC
  double f1 = eop[day_before][5];      // first value of UT1-UTC
  double f2 = eop[day_before + 1][5];  // second value of UT1-UTC
  double ut1 = (utc - day_before) * f2 - (utc - day_before - 1) * f1;
  return ut1;
}

double TT2UTC(double mjd_tt) {
  double sec_tt_utc = 32.184 + 37.0;  // valid from 01/01/2017 onwards.
  double mdj_utc = mjd_tt - sec_tt_utc / 86400.0;
  return mdj_utc;
}

double MeanObliquity(double mjd_TT) {
  const double T = (mjd_TT - MJD_J2000) / 36525.0;

  return DEG_TO_RAD * (23.43929111 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T / 3600.0);
}

double EqEquinoxes(double mjd_TT) {
  double dpsi, deps;  // Nutation angles

  NutationAngles(mjd_TT, dpsi, deps);
  static int count = 0;
  if (count == 0) {
    count++;
    std::cout << "dpsi = " << dpsi << std::endl;
    std::cout << "meanObli = " << MeanObliquity(mjd_TT) << std::endl;
  }
  // Equation of the equinoxes
  return dpsi * cos(MeanObliquity(mjd_TT));
};

double GMST(double mjd_UT1) {
  // Variables
  double mjd_UT1_0, UT1, T_0, T, gmst;

  // Mean Sidereal Time
  mjd_UT1_0 = floor(mjd_UT1);
  static int count = 0;
  UT1 = 86400.0 * (mjd_UT1 - mjd_UT1_0);  // seconds in UT1 of that day [s]
  T_0 = (mjd_UT1_0 - MJD_J2000) / 36525.0;
  T = (mjd_UT1 - MJD_J2000) / 36525.0;

  // gmst = 24110.54841 + 8640184.812866 * T + 1.002737909350795 * UT1 + (0.093104 - 6.2e-6 * T) * T * T;  // [s] - How I think it should be (https://gssc.esa.int/navipedia/index.php/CEP_to_ITRF)
  gmst = 24110.54841 + 8640184.812866 * T_0 + 1.002737909350795 * UT1 + (0.093104 - 6.2e-6 * T) * T * T;  // [s] - Montenbruck

  if (count == 0) {
    count++;
    // change precision output
    std::cout.precision(15);
    std::cout << "mjd_UT1 = " << mjd_UT1 << std::endl;
    std::cout << "mjd_UT1_0 = " << mjd_UT1_0 << std::endl;
    std::cout << "mjd_UT1 - mjd_UT1_0 = " << mjd_UT1 - mjd_UT1_0 << std::endl;
    std::cout << "UT1 = " << UT1 << std::endl;
    std::cout << "T_0 = " << T_0 << std::endl;
    std::cout << "T = " << T << std::endl;
    std::cout << "gmst = " << gmst << std::endl;
    std::cout << "gmst_alternative = " << modf(gmst / 86400.0, &gmst) * 2 * PI << std::endl;
  }
  // first we extract the fractional part of the gmst (which is the fraction of a day (or revolution)), and then we convert it to radians. The passing reference gmst is just to avoid a warning, but here we don't need it (it would be the integer part of gmst).
  return modf(gmst / 86400.0, &gmst) * 2 * PI;  // [rad] - Correction based on Montenbruck
  // return fmod(gmst * SEC_TO_RAD, 2 * PI);
}

double GAST(double mjd_TT) {
  // we assume UT1 = UTC
  double incr = (32.184 + 37.0) / 86400.0;  // TT-UT1 [days]
  double mjd_UT1 = mjd_TT - incr;
  static int count = 0;
  if (count == 0) {
    count++;
    std::cout << "equinox = " << EqEquinoxes(mjd_TT) << std::endl;
    std::cout << "gmstGAST = " << GMST(mjd_TT) << std::endl;
    std::cout << "mjd_TTGAST = " << mjd_TT << std::endl;
  }
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
  return floor(365.25 * (yyyy + 4716)) + floor(30.6001 * (mm + 1)) + dd + b - 1524.5 - 2400000.5;
}

Matrix PrecessionMatrix(double mjd_TT) {
  // Constants
  const double T = (mjd_TT - MJD_J2000) / 36525.0;

  // Variables
  double zeta, z, theta;

  // Precession angles
  zeta = (2306.2181 + (0.30188 + 0.017998 * T) * T) * T * SEC_TO_RAD;
  z = zeta + (0.79280 + 0.000205 * T) * T * T * SEC_TO_RAD;
  theta = (2004.3109 - (0.42665 + 0.041833 * T) * T) * T * SEC_TO_RAD;

  static int count = 0;

  if (count == 0) {
    count++;
    std::cout << "MJD_TT: " << mjd_TT << std::endl;
    std::cout << "T: " << T << std::endl;
    std::cout << "zeta: " << zeta << std::endl;
    std::cout << "z: " << z << std::endl;
    std::cout << "theta: " << theta << std::endl;
  }

  // Precession matrix
  Matrix R(3, -z), S(2, theta), V(3, -zeta);
  return R * S * V;
}

void NutationAngles(double mjd_TT, double& dpsi, double& deps) {
  // Constants
  const double T = (mjd_TT - MJD_J2000) / 36525.0;
  const double T2 = T * T;
  const double T3 = T2 * T;
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

  l = fmod(485866.733 + (1325.0 * rev + 715922.633) * T + 31.310 * T2 + 0.064 * T3, rev);
  lp = fmod(1287099.804 + (99.0 * rev + 1292581.224) * T - 0.577 * T2 - 0.012 * T3, rev);
  F = fmod(335778.877 + (1342.0 * rev + 295263.137) * T - 13.257 * T2 + 0.011 * T3, rev);
  D = fmod(1072261.307 + (1236.0 * rev + 1105601.328) * T - 6.891 * T2 + 0.019 * T3, rev);
  Om = fmod(450160.280 - (5.0 * rev + 482890.539) * T + 7.455 * T2 + 0.008 * T3, rev);

  // Nutation in longitude and obliquity [rad]

  deps = dpsi = 0.0;
  for (int i = 0; i < nut_coeffs; i++) {
    arg = (C[i][0] * l + C[i][1] * lp + C[i][2] * F + C[i][3] * D + C[i][4] * Om) * SEC_TO_RAD;
    dpsi += (C[i][5] + C[i][6] * T) * sin(arg);
    deps += (C[i][7] + C[i][8] * T) * cos(arg);
  };
  dpsi = 1.0E-5 * dpsi * SEC_TO_RAD;
  deps = 1.0E-5 * deps * SEC_TO_RAD;
}

Matrix NutationMatrix(double mjd_TT) {
  double dpsi, deps, eps;

  // Mean obliquity of the ecliptic
  eps = MeanObliquity(mjd_TT);

  // Nutation in longitude and obliquity
  NutationAngles(mjd_TT, dpsi, deps);

  // Transformation from mean to true equator and equinox
  Matrix R(1, -eps - deps), S(3, dpsi), T(1, eps);
  return R * S * T;
}

Matrix RotationMatrix(double mjd_TT) {
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

  xp *= SEC_TO_RAD;
  yp *= SEC_TO_RAD;
}

Matrix PolarMotionMatrix(double mjd_TT) {
  double xp, yp;

  // Polar motion
  PolarMotion(mjd_TT, xp, yp);

  // Transformation from ITRS to GCRS
  Matrix R(2, -xp), S(1, -yp);
  return R * S;
}

Vector Sun(double Mjd_TT) {
  // Constants

  const double eps = 23.43929111 * DEG_TO_RAD;      // Obliquity of J2000 ecliptic
  const double T = (Mjd_TT - MJD_J2000) / 36525.0;  // Julian cent. since J2000

  // Variables

  double L, M, r;
  Vector r_Sun(3);

  // // Mean anomaly, ecliptic longitude and radius

  // M = pi2 * Frac(0.9931267 + 99.9973583 * T);  // [rad]
  // L = pi2 * Frac(0.7859444 + M / pi2 + (6892.0 * sin(M) + 72.0 * sin(2.0 * M)) / 1296.0e3);  // [rad]
  // r = 149.619e9 - 2.499e9 * cos(M) - 0.021e9 * cos(2 * M);             // [m]

  // // Equatorial position vector

  // r_Sun = R_x(-eps) * Vector(r * cos(L), r * sin(L), 0.0);

  return r_Sun;
}

Vector Moon(double Mjd_TT) {
  // Constants

  const double eps = 23.43929111 * DEG_TO_RAD;      // Obliquity of J2000 ecliptic
  const double T = (Mjd_TT - MJD_J2000) / 36525.0;  // Julian cent. since J2000

  // Variables

  double L_0, l, lp, F, D, dL, S, h, N;
  double L, B, R, cosB;
  Vector r_Moon(3);

  // // Mean elements of lunar orbit

  // L_0 = Frac(0.606433 + 1336.851344 * T);      // Mean longitude [rev]
  //                                              // w.r.t. J2000 equinox
  // l = pi2 * Frac(0.374897 + 1325.552410 * T);  // Moon's mean anomaly [rad]
  // lp = pi2 * Frac(0.993133 + 99.997361 * T);   // Sun's mean anomaly [rad]
  // D = pi2 * Frac(0.827361 + 1236.853086 * T);  // Diff. long. Moon-Sun [rad]
  // F = pi2 * Frac(0.259086 + 1342.227825 * T);  // Argument of latitude

  // // Ecliptic longitude (w.r.t. equinox of J2000)

  // dL = +22640 * sin(l) - 4586 * sin(l - 2 * D) + 2370 * sin(2 * D) + 769 * sin(2 * l) - 668 * sin(lp) - 412 * sin(2 * F) - 212 * sin(2 * l - 2 * D) - 206 * sin(l + lp - 2 * D) + 192 * sin(l + 2 * D) - 165 * sin(lp - 2 * D) - 125 * sin(D) - 110 * sin(l + lp) + 148 * sin(l - lp) - 55 * sin(2 * F - 2 * D);

  // L = pi2 * Frac(L_0 + dL / 1296.0e3);  // [rad]

  // // Ecliptic latitude

  // S = F + (dL + 412 * sin(2 * F) + 541 * sin(lp)) / Arcs;
  // h = F - 2 * D;
  // N = -526 * sin(h) + 44 * sin(l + h) - 31 * sin(-l + h) - 23 * sin(lp + h) + 11 * sin(-lp + h) - 25 * sin(-2 * l + F) + 21 * sin(-l + F);

  // B = (18520.0 * sin(S) + N) / Arcs;  // [rad]

  // cosB = cos(B);

  // // Distance [m]

  // R = 385000e3 - 20905e3 * cos(l) - 3699e3 * cos(2 * D - l) - 2956e3 * cos(2 * D) - 570e3 * cos(2 * l) + 246e3 * cos(2 * l - 2 * D) - 205e3 * cos(lp - 2 * D) - 171e3 * cos(l + 2 * D) - 152e3 * cos(l + lp - 2 * D);

  // // Equatorial coordinates

  // r_Moon = R_x(-eps) * Vector(R * cos(L) * cosB, R * sin(L) * cosB, R * sin(B));

  return r_Moon;
}
