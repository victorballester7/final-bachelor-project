#include "../include/ReferenceSystems.h"

#include "../include/LinearAlgebra.h"
#include "../include/constants.h"

Matrix GCRF2Perifocal(double i, double Omega, double omega) {
  Matrix R(3, Omega), S(1, i), T(3, omega);
  return T * S * R;
}

Matrix GCRF2ITRF(double mjd_TT) {
  double mjd_UT1 = mjd_TT;  // we assume that UT1 = TT
  return PolarMotionMatrix(mjd_TT) * RotationMatrix(mjd_UT1) * NutationMatrix(mjd_TT) * PrecessionMatrix(mjd_TT);
}

// double TT2JDcenturies_TT(double t) {
//   // TT (Terrestrial Time) Julian Date
//   return Mjd_TT + 2400000.5;
// }

double MeanObliquity(double mjd_TT) {
  const double T = (mjd_TT - MJD_J2000) / 36525.0;

  return DEG_TO_RAD * (23.43929111 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T / 3600.0);
}

double EqEquinoxes(double mjd_TT) {
  double dpsi, deps;  // Nutation angles

  NutationAngles(mjd_TT, dpsi, deps);

  // Equation of the equinoxes
  return dpsi * cos(MeanObliquity(mjd_TT));
};

double GMST(double mjd_UT1) {
  // Variables
  double mjd_0, UT1, T_0, T, gmst;

  // Mean Sidereal Time
  mjd_0 = floor(mjd_UT1);
  UT1 = 86400.0 * (mjd_UT1 - mjd_0);  // seconds in UT1 of that day [s]
  T_0 = (mjd_0 - MJD_J2000) / 36525.0;
  T = (mjd_UT1 - MJD_J2000) / 36525.0;

  gmst = 24110.54841 + 8640184.812866 * T_0 + 1.002737909350795 * UT1 + (0.093104 - 6.2e-6 * T) * T * T;  // [s]

  return fmod(gmst * SEC_TO_RAD, 2 * PI);
}

double GAST(double mjd_UT1) {
  // we asume that mjd_UT1 = mjd_TT when passed to the EqEquinoxes function
  return fmod(GMST(mjd_UT1) + EqEquinoxes(mjd_UT1), 2 * PI);
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

double get_mjd_TT(double yyyy, double mm, double dd) {
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

  const int N_coeff = 106;
  const long C[N_coeff][9] =
      {
          //
          // l  l' F  D Om    dpsi    *T     deps     *T       #
          //
          {0, 0, 0, 0, 1, -1719960, -1742, 920250, 89},  //   1
          {0, 0, 0, 0, 2, 20620, 2, -8950, 5},           //   2
          {-2, 0, 2, 0, 1, 460, 0, -240, 0},             //   3
          {2, 0, -2, 0, 0, 110, 0, 0, 0},                //   4
          {-2, 0, 2, 0, 2, -30, 0, 10, 0},               //   5
          {1, -1, 0, -1, 0, -30, 0, 0, 0},               //   6
          {0, -2, 2, -2, 1, -20, 0, 10, 0},              //   7
          {2, 0, -2, 0, 1, 10, 0, 0, 0},                 //   8
          {0, 0, 2, -2, 2, -131870, -16, 57360, -31},    //   9
          {0, 1, 0, 0, 0, 14260, -34, 540, -1},          //  10
          {0, 1, 2, -2, 2, -5170, 12, 2240, -6},         //  11
          {0, -1, 2, -2, 2, 2170, -5, -950, 3},          //  12
          {0, 0, 2, -2, 1, 1290, 1, -700, 0},            //  13
          {2, 0, 0, -2, 0, 480, 0, 10, 0},               //  14
          {0, 0, 2, -2, 0, -220, 0, 0, 0},               //  15
          {0, 2, 0, 0, 0, 170, -1, 0, 0},                //  16
          {0, 1, 0, 0, 1, -150, 0, 90, 0},               //  17
          {0, 2, 2, -2, 2, -160, 1, 70, 0},              //  18
          {0, -1, 0, 0, 1, -120, 0, 60, 0},              //  19
          {-2, 0, 0, 2, 1, -60, 0, 30, 0},               //  20
          {0, -1, 2, -2, 1, -50, 0, 30, 0},              //  21
          {2, 0, 0, -2, 1, 40, 0, -20, 0},               //  22
          {0, 1, 2, -2, 1, 40, 0, -20, 0},               //  23
          {1, 0, 0, -1, 0, -40, 0, 0, 0},                //  24
          {2, 1, 0, -2, 0, 10, 0, 0, 0},                 //  25
          {0, 0, -2, 2, 1, 10, 0, 0, 0},                 //  26
          {0, 1, -2, 2, 0, -10, 0, 0, 0},                //  27
          {0, 1, 0, 0, 2, 10, 0, 0, 0},                  //  28
          {-1, 0, 0, 1, 1, 10, 0, 0, 0},                 //  29
          {0, 1, 2, -2, 0, -10, 0, 0, 0},                //  30
          {0, 0, 2, 0, 2, -22740, -2, 9770, -5},         //  31
          {1, 0, 0, 0, 0, 7120, 1, -70, 0},              //  32
          {0, 0, 2, 0, 1, -3860, -4, 2000, 0},           //  33
          {1, 0, 2, 0, 2, -3010, 0, 1290, -1},           //  34
          {1, 0, 0, -2, 0, -1580, 0, -10, 0},            //  35
          {-1, 0, 2, 0, 2, 1230, 0, -530, 0},            //  36
          {0, 0, 0, 2, 0, 630, 0, -20, 0},               //  37
          {1, 0, 0, 0, 1, 630, 1, -330, 0},              //  38
          {-1, 0, 0, 0, 1, -580, -1, 320, 0},            //  39
          {-1, 0, 2, 2, 2, -590, 0, 260, 0},             //  40
          {1, 0, 2, 0, 1, -510, 0, 270, 0},              //  41
          {0, 0, 2, 2, 2, -380, 0, 160, 0},              //  42
          {2, 0, 0, 0, 0, 290, 0, -10, 0},               //  43
          {1, 0, 2, -2, 2, 290, 0, -120, 0},             //  44
          {2, 0, 2, 0, 2, -310, 0, 130, 0},              //  45
          {0, 0, 2, 0, 0, 260, 0, -10, 0},               //  46
          {-1, 0, 2, 0, 1, 210, 0, -100, 0},             //  47
          {-1, 0, 0, 2, 1, 160, 0, -80, 0},              //  48
          {1, 0, 0, -2, 1, -130, 0, 70, 0},              //  49
          {-1, 0, 2, 2, 1, -100, 0, 50, 0},              //  50
          {1, 1, 0, -2, 0, -70, 0, 0, 0},                //  51
          {0, 1, 2, 0, 2, 70, 0, -30, 0},                //  52
          {0, -1, 2, 0, 2, -70, 0, 30, 0},               //  53
          {1, 0, 2, 2, 2, -80, 0, 30, 0},                //  54
          {1, 0, 0, 2, 0, 60, 0, 0, 0},                  //  55
          {2, 0, 2, -2, 2, 60, 0, -30, 0},               //  56
          {0, 0, 0, 2, 1, -60, 0, 30, 0},                //  57
          {0, 0, 2, 2, 1, -70, 0, 30, 0},                //  58
          {1, 0, 2, -2, 1, 60, 0, -30, 0},               //  59
          {0, 0, 0, -2, 1, -50, 0, 30, 0},               //  60
          {1, -1, 0, 0, 0, 50, 0, 0, 0},                 //  61
          {2, 0, 2, 0, 1, -50, 0, 30, 0},                //  62
          {0, 1, 0, -2, 0, -40, 0, 0, 0},                //  63
          {1, 0, -2, 0, 0, 40, 0, 0, 0},                 //  64
          {0, 0, 0, 1, 0, -40, 0, 0, 0},                 //  65
          {1, 1, 0, 0, 0, -30, 0, 0, 0},                 //  66
          {1, 0, 2, 0, 0, 30, 0, 0, 0},                  //  67
          {1, -1, 2, 0, 2, -30, 0, 10, 0},               //  68
          {-1, -1, 2, 2, 2, -30, 0, 10, 0},              //  69
          {-2, 0, 0, 0, 1, -20, 0, 10, 0},               //  70
          {3, 0, 2, 0, 2, -30, 0, 10, 0},                //  71
          {0, -1, 2, 2, 2, -30, 0, 10, 0},               //  72
          {1, 1, 2, 0, 2, 20, 0, -10, 0},                //  73
          {-1, 0, 2, -2, 1, -20, 0, 10, 0},              //  74
          {2, 0, 0, 0, 1, 20, 0, -10, 0},                //  75
          {1, 0, 0, 0, 2, -20, 0, 10, 0},                //  76
          {3, 0, 0, 0, 0, 20, 0, 0, 0},                  //  77
          {0, 0, 2, 1, 2, 20, 0, -10, 0},                //  78
          {-1, 0, 0, 0, 2, 10, 0, -10, 0},               //  79
          {1, 0, 0, -4, 0, -10, 0, 0, 0},                //  80
          {-2, 0, 2, 2, 2, 10, 0, -10, 0},               //  81
          {-1, 0, 2, 4, 2, -20, 0, 10, 0},               //  82
          {2, 0, 0, -4, 0, -10, 0, 0, 0},                //  83
          {1, 1, 2, -2, 2, 10, 0, -10, 0},               //  84
          {1, 0, 2, 2, 1, -10, 0, 10, 0},                //  85
          {-2, 0, 2, 4, 2, -10, 0, 10, 0},               //  86
          {-1, 0, 4, 0, 2, 10, 0, 0, 0},                 //  87
          {1, -1, 0, -2, 0, 10, 0, 0, 0},                //  88
          {2, 0, 2, -2, 1, 10, 0, -10, 0},               //  89
          {2, 0, 2, 2, 2, -10, 0, 0, 0},                 //  90
          {1, 0, 0, 2, 1, -10, 0, 0, 0},                 //  91
          {0, 0, 4, -2, 2, 10, 0, 0, 0},                 //  92
          {3, 0, 2, -2, 2, 10, 0, 0, 0},                 //  93
          {1, 0, 2, -2, 0, -10, 0, 0, 0},                //  94
          {0, 1, 2, 0, 1, 10, 0, 0, 0},                  //  95
          {-1, -1, 0, 2, 1, 10, 0, 0, 0},                //  96
          {0, 0, -2, 0, 1, -10, 0, 0, 0},                //  97
          {0, 0, 2, -1, 2, -10, 0, 0, 0},                //  98
          {0, 1, 0, 2, 0, -10, 0, 0, 0},                 //  99
          {1, 0, -2, -2, 0, -10, 0, 0, 0},               // 100
          {0, -1, 2, 0, 1, -10, 0, 0, 0},                // 101
          {1, 1, 0, -2, 1, -10, 0, 0, 0},                // 102
          {1, 0, -2, 2, 0, -10, 0, 0, 0},                // 103
          {2, 0, 0, 2, 0, 10, 0, 0, 0},                  // 104
          {0, 0, 2, 4, 2, -10, 0, 0, 0},                 // 105
          {0, 1, 0, 1, 0, 10, 0, 0, 0}                   // 106
      };

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
  for (int i = 0; i < N_coeff; i++) {
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

Matrix RotationMatrix(double mjd_UT1) {
  return Matrix(3, GAST(mjd_UT1));
}

void PolarMotion(double mjd_TT, double& xp, double& yp) {
  xp = 0 * mjd_TT;
  yp = 0;
}

Matrix PolarMotionMatrix(double mjd_TT) {
  double xp, yp;

  // Polar motion
  PolarMotion(mjd_TT, xp, yp);

  // Transformation from ITRS to GCRS
  Matrix R(2, -xp), S(1, -yp);
  return R * S;
}