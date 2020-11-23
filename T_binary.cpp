#include "Triangular.h"

int main() {
  START();

  int M1 = 9;
  int M2 = M1 + 1;
  double Js[3] = {3.0, 1.0, 1.0};
  double EPS = 1e-12;

  Triangular T1(Js, M1, EPS);
  double* vo_1 = alloc_dvector(T1.dim);
  double* vn_1 = alloc_dvector(T1.dim);
  double* vtmp1_1 = alloc_dvector(T1.dim2);
  double* vtmp2_1 = alloc_dvector(T1.dim2);
  double* v1R_1 = alloc_dvector(T1.dim);
  double* v1L_1 = alloc_dvector(T1.dim);
  double* v2R_1 = alloc_dvector(T1.dim);

  Triangular T2(Js, M2, EPS);
  double* vo_2 = alloc_dvector(T2.dim);
  double* vn_2 = alloc_dvector(T2.dim);
  double* vtmp1_2 = alloc_dvector(T2.dim2);
  double* vtmp2_2 = alloc_dvector(T2.dim2);
  double* v1R_2 = alloc_dvector(T2.dim);
  double* v1L_2 = alloc_dvector(T2.dim);
  double* v2R_2 = alloc_dvector(T2.dim);

  double l = 0.2, r = 0.7, c;
  while (r - l > EPS) {
    c = (l + r) / 2.0;
    double xi1 =
        T1.calc_xi(c, vo_1, vn_1, vtmp1_1, vtmp2_1, v1R_1, v1L_1, v2R_1);
    double xi2 =
        T2.calc_xi(c, vo_2, vn_2, vtmp1_2, vtmp2_2, v1R_2, v1L_2, v2R_2);
    (M1 / xi1 > M2 / xi2 ? l : r) = c;
  }
  double Tc = (l + r) / 2;
  printf("M = %d, %d\nTc = %.11f\n", M1, M2, Tc);
  printf("Tc (exact) = %.11f\n", T1.exact_Tc());

  END();
}