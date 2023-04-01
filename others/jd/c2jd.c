#include <stdio.h>
#include "jd.h"

int main (int argc, char *argv[]) {
   double d, m, y;
   if (argc<4
	 || sscanf(argv[1], "%lf", &y)!=1
	 || sscanf(argv[2], "%lf", &m)!=1
	 || sscanf(argv[3], "%lf", &d)!=1
	 ) {
      fprintf(stderr, "%s y m d\n", argv[0]);
      return -1;
   }
   printf("%.16G\n", c2jd(y,m,d));
   return 0;
}
