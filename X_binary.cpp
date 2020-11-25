#include "Xsquare.h"

int main() {
  START();

  /*
    J = (1,1,1,1)
    M = 13,14
    Tc = 0.65704177490
  */

  int M1 = 4, M2 = M1 + 1;
  double Js[4] = {1.0, 1.0, 1.0, 1.0};
  double EPS = 1e-12;
  double L = 0.4, R = 0.67;

  Xsquare X1(Js, M1, EPS);
  double* vo_1 = alloc_dvector(X1.dim);
  double* vn_1 = alloc_dvector(X1.dim);
  double* vtmp1_1 = alloc_dvector(X1.dimt);
  double* vtmp2_1 = alloc_dvector(X1.dimt);
  double* v1R_1 = alloc_dvector(X1.dim);
  double* v1L_1 = alloc_dvector(X1.dim);
  double* v2R_1 = alloc_dvector(X1.dim);

  Xsquare X2(Js, M2, EPS);
  double* vo_2 = alloc_dvector(X2.dim);
  double* vn_2 = alloc_dvector(X2.dim);
  double* vtmp1_2 = alloc_dvector(X2.dimt);
  double* vtmp2_2 = alloc_dvector(X2.dimt);
  double* v1R_2 = alloc_dvector(X2.dim);
  double* v1L_2 = alloc_dvector(X2.dim);
  double* v2R_2 = alloc_dvector(X2.dim);

  double xi1, xi2;
  double l = L, r = R, c;
  while (r - l > EPS) {
    c = (l + r) / 2.0;
    xi1 = X1.calc_xi(c, vo_1, vn_1, vtmp1_1, vtmp2_1, v1R_1, v1L_1, v2R_1);
    xi2 = X2.calc_xi(c, vo_2, vn_2, vtmp1_2, vtmp2_2, v1R_2, v1L_2, v2R_2);

    if (xi1 == -1.0 || xi2 == -1.0)
      r = c;
    else
      ((M1 / xi1 > M2 / xi2) ? l : r) = c;
  }

  double Tc = (l + r) / 2.0;
  if (xi1 == -1.0 || xi2 == -1.0)
    Tc = -1.0;

  printf("J = (%.1f, %.1f, %.1f, %.1f)\n", Js[0], Js[1], Js[2], Js[3]);
  printf("M = %d, %d\n%.11f\n", M1, M2, Tc);

  if (Tc != -1.0) {
    // increase
    double next_L = c, next_R = 2 * c - L;
    printf("Inc\ndouble L = %.11f, R = %.11f;\n", next_L, next_R);

    // decrease
    next_R = c, next_L = 2 * c - R;
    printf("Dec\ndouble L = %.11f, R = %.11f;\n", next_L, next_R);
  }

  END();
}