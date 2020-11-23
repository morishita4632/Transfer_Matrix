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
  double* vtmp1 = alloc_dvector(T.dim2);
  double* vtmp2 = alloc_dvector(T.dim2);
  double* v1R = alloc_dvector(T.dim);
  double* v1L = alloc_dvector(T.dim);
  double* v2R = alloc_dvector(T.dim);

  double xi = T.calc_xi(temperature, vo, vn, vtmp1, vtmp2, v1R, v1L, v2R);
  if (xi == -1.0)
    printf("NOT ORTHOGONAL!\n");
  else
    printf("xi = %.9f\n", xi);

  END();
}