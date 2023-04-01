#define MJD2000_ORIG 2451544.5

/*
 * Altres origens habituals, a confirmar:
 * - 17.0-Nov-1858 : JD 2400000.5
 * - 1.0-Gen-1950  : JD 2433282.5
 */

/* Canvi de data de calendari a data juliana */
double c2jd (double y, double m, double d);

/* Canvi de data juliana a data de calendari */
void jd2c (double jd, double *yy, double *mm, double *dd);
