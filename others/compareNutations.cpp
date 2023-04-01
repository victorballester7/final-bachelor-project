
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <iostream>

#include "TT_to_TDB.h"
#include "jpleph/jpleph.h"
#include "monte/SAT_Const.h"
#include "monte/SAT_RefSys.h"

int main(int argc, char const *argv[]) {
  double jd_TT_large, jd_TT_small;
  double mjd_TT;      // modified Julian Date in TT
  double dpsi, deps;  // longitude and latitude (obliquity) of the nutation
  struct jpleph jph;

  if (argc != 3 || sscanf(argv[1], "%lf", &jd_TT_large) != 1 || sscanf(argv[2], "%lf", &jd_TT_small) != 1) {
    printf("Execute as ./compareNutations JD_TT_large JD_TT_small\n");
    return -1;
  }
  mjd_TT = (jd_TT_large - 2400000.5) + jd_TT_small;
  NutAngles(mjd_TT, dpsi, deps);

  printf("\nTerrestrial Time = %lf\n\n\n", jd_TT_large + jd_TT_small);
  printf("Nutation angles Montenbruck:\n\tdpsi = %.16lf  [arcsec]\n\tdeps = %.16lf  [arcsec]\n\n", dpsi * Arcs, deps * Arcs);

  char filename[] = "jpleph/de405.lnx";
  assert(jpleph_open(filename, &jph, 0 /*rflag*/) >= 0);
  double d0, dlen;
  int ncfs;
  double cfs[jph.mxlncfs];
  jpleph_cfs(&jph, jd_TT_large + jd_TT_small /*d*/, 12, 1, &d0, &dlen, &ncfs, cfs);
  dpsi = jpleph_cfs_aval(jd_TT_large, jd_TT_small, d0, dlen, ncfs, cfs);
  jpleph_cfs(&jph, jd_TT_large + jd_TT_small /*d*/, 12, 2, &d0, &dlen, &ncfs, cfs);
  deps = jpleph_cfs_aval(jd_TT_large, jd_TT_small, d0, dlen, ncfs, cfs);
  printf("Nutation angles JPL Ephemerids (using JD_TT, which is incorrect):\n\tdpsi = %.16lf  [arcsec]\n\tdeps = %.16lf  [arcsec]\n\n", dpsi * Arcs, deps * Arcs);

  double jd_TDB_large = jd_TT_large, jd_TDB_small = jd_TT_small;
  jd_TDB_small += TT_to_TDB(jd_TT_large, jd_TT_small);
  jpleph_cfs(&jph, jd_TDB_large + jd_TDB_small /*d*/, 12, 1, &d0, &dlen, &ncfs, cfs);
  dpsi = jpleph_cfs_aval(jd_TDB_large, jd_TDB_small, d0, dlen, ncfs, cfs);
  jpleph_cfs(&jph, jd_TDB_large + jd_TDB_small /*d*/, 12, 2, &d0, &dlen, &ncfs, cfs);
  deps = jpleph_cfs_aval(jd_TDB_large, jd_TDB_small, d0, dlen, ncfs, cfs);
  printf("Nutation angles JPL Ephemerids (using JD_TDB):\n\tdpsi = %.16lf  [arcsec]\n\tdeps = %.16lf  [arcsec]\n\n", dpsi * Arcs, deps * Arcs);
  return 0;
}
