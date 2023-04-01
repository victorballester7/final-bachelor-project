#include <stdio.h>
#include "jd.h"

int main (int argc, char *argv[]) {
   double jd, y, m, d;
   if (argc!=2
	 || sscanf(argv[1], "%lf", &jd)!=1
	 ) {
      fprintf(stderr, "%s jd\n", argv[0]);
      return -1;
   }
   jd2c(jd,&y,&m,&d);
   printf("%.16G %.16G %.16G\n", y, m, d);
   return 0;
}
