#include <cassert>
#include <cstdio>

#include "SAT_Force.h"
#include "SAT_VecMat.h"

#define M 3
#define N 6

int main(void) {
  Vector vfM(-738.672, 4500.57, 6064.25), rM(6.83509e+06, 1.10735e+06, 9550.44);
  double id[M * M] = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
  Matrix E(id, M, M);
  /* Cridem Monte */
  // vfM = AccelHarmonic(rM, E, Grav.GM, Grav.R_ref, Grav.CS, 8 /*n_max*/, 8 /*m_max*/);
  double mjd_tt = 59945.33310694993;
  int n_max = 8;
  int m_max = n_max;
  vfM = AccelMainCustom(mjd_tt, rM, vfM, n_max, m_max, 0, 0, 0, 0);
  printf("Monte tallat a %dx%d:\n\t%.16G %.16G %.16G\n", n_max, m_max, vfM(0), vfM(1), vfM(2));
  return 0;
}
