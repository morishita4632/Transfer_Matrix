#include "Xsquare.h"

int main() {
  START();

  int M = 16;
  double Js[4] = {1.0, 1.0, 1.0, 1.0};
  double temperature = 1.0;
  double EPS = 1e-12;

  Xsquare X(Js, M, EPS);
  double* vo = alloc_dvector(X.dim);
  double* vn = alloc_dvector(X.dim);
  double* vtmp1 = alloc_dvector(X.dimt);
  double* vtmp2 = alloc_dvector(X.dimt);
  double* v1R = alloc_dvector(X.dim);
  double* v1L = alloc_dvector(X.dim);
  double* v2R = alloc_dvector(X.dim);

  double lmd1 = X.power1_R(temperature, vo, vn, v1R, vtmp1, vtmp2);
  X.power1_L(temperature, vo, vn, v1L, vtmp1, vtmp2);
  double lmd2 =
      X.power2(temperature, lmd1, vo, vn, v1R, v1L, vtmp1, vtmp2, v2R);
  double xi = X.calc_xi(temperature, vo, vn, vtmp1, vtmp2, v1R, v1L, v2R);

  printf("λ1 = %.12f\n", lmd1);
  printf("λ2 = %.12f\n", lmd2);
  printf("ξ = %.12f\n", xi);

  END();
}