#include <cassert>
#include <cstdio>

#include "SAT_Force.h"
#include "SAT_VecMat.h"

#define M 3
#define N 6

int main(void) {
  Vector vfM(M), rM(1.e7, 1.e7, 1.e7);
  double id[M * M] = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
  Matrix E(id, M, M);
  /* Cridem Monte */
  vfM = AccelHarmonic(rM, E, Grav.GM, Grav.R_ref, Grav.CS, 8 /*n_max*/, 8 /*m_max*/);
  printf("Monte tallat a 8x8:\n\t%.16G %.16G %.16G\n",
         vfM(0), vfM(1), vfM(2));
  return 0;
}
