#include "Triangular.h"

int main() {
  START();

  int M = 6;
  double Js[3] = {7.0, 1.0, 1.0};
  double EPS = 1e-12;
  double temperature = 0.6;

  Triangular T(Js, M, EPS);
  double* vo = alloc_dvector(T.dim);
  double* vn = alloc_dvector(T.dim);
  double* vtmp1 = alloc_dvector(T.dimt);
  double* vtmp2 = alloc_dvector(T.dimt);
  double* v1 = alloc_dvector(T.dim);
  double* v2 = alloc_dvector(T.dim);

  double lmd1 = T.power1(temperature, vo, vn, v1, vtmp1, vtmp2);
  double lmd2 = T.power2(temperature, lmd1, vo, vn, v1, v2, vtmp1, vtmp2);
  double xi = T.calc_xi(temperature, vo, vn, vtmp1, vtmp2, v1, v2);

  printf("ξ = %.12f\n", xi);
  printf("λ1 = %.12f\n", lmd1);
  printf("λ2 = %.12f\n", lmd2);

  END();
}