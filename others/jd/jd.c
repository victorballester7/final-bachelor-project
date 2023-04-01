#include <stdio.h>
#include <math.h>

/* Canvi de data de calendari a data juliana */
double c2jd (double y, double m, double d) {
   double a, b;
   if (m<=2) { y--; m+=12; }
   a=floor(y/100);
   b=y+0.01*m+0.0001*d>=1582.1015 ?
      2-a+floor(a/4)	/* correccio per la reforma gregoriana */
      : 0;
   return floor(365.25*(y+4716))+floor(30.6001*(m+1))+d+b-1524.5;
}

/* Canvi de data juliana a data de calendari */
void jd2c (double jd, double *yy, double *mm, double *dd) {
   double z, f, b, c, d, e;
   jd+=0.5;
   z=floor(jd);
   f=jd-z;
   b=z;
   if (z>=2299161) {	/* correccio per la reforma gregoriana */
      double alf=floor((z-1867126.25)/36524.25);
      b+=1+alf-floor(alf/4);
   }
   b+=1524;
   c=floor((b-122.1)/365.25);
   d=floor(365.25*c);
   e=floor((b-d)/30.6001);
   *dd=b-d-floor(30.6001*e)+f;
   *mm=e<14 ? e-1 : e-13;
   *yy=4716;
   if (*mm<=2) (*yy)--;
   *yy=c-(*yy);
}
