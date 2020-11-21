#include "X_Square.h"

int main() {
  int M = 16;
  double Js[4] = {1.0, 1.0, 1.0, 0.0};
  double temperature = 1.0;
  double EPS = 1e-12;

  X_Square X(Js, M, EPS);
  double* vo = alloc_dvector(X.dim);
  double* vn = alloc_dvector(X.dim);
  double* vtmp1 = alloc_dvector(X.dim4);
  double* vtmp2 = alloc_dvector(X.dim4);
  double* v1R = alloc_dvector(X.dim);
  double* v1L = alloc_dvector(X.dim);


  double lmd1 = X.power1_R(temperature, vo, vn, v1R, vtmp1, vtmp2);
  printf("%.9f\n", lmd1);
  // X.power1_L(temperature, vo, vn, v1L, vtmp1, vtmp2);
  // double lmd2 = X.power2(temperature, lmd1, vo, vn, v1R, v1L, vtmp1, vtmp2);
  // double xi = X.calc_xi(temperature, vo, vn, vtmp1, vtmp2, v1R, v1L);

  // printf("λ1 = %.12f\n", lmd1);
  // printf("λ2 = %.12f\n", lmd2);
  // printf("ξ = %.12f\n", xi);
}