#include "Triangular.h"

// v2Rにはもともと0.0の成分があるのでそこを無視する必要がある

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
  double* v1R = alloc_dvector(T.dim);
  double* v1L = alloc_dvector(T.dim);
  double* v2R = alloc_dvector(T.dim);

  double lmd1 = T.power1_R(temperature, vo, vn, v1R, vtmp1, vtmp2);
  T.power1_L(temperature, vo, vn, v1L, vtmp1, vtmp2);
  double inner = dot(v1R, v1L, T.dim);
  for (int i = 0; i < T.dim; i++)
    v1L[i] /= inner;
  double lmd2 =
      T.power2(temperature, lmd1, vo, vn, v1R, v1L, v2R, vtmp1, vtmp2);

  int cnt;
  double mult;
  double test_eps = 1e-5;

  printf("λ1 = %.12f\n", lmd1);
  vcopy(vtmp1, v1R, T.dim);
  T.product(temperature, vtmp1, vtmp2);
  cnt = 0;
  for (int i = 0; i < T.dim; i++) {
    mult = vtmp1[i] / v1R[i];
    cnt += abs(mult - lmd1) > test_eps;
  }
  printf("%d\n", cnt);

  vcopy(vtmp1, v1L, T.dim);
  T.product(temperature, vtmp1, vtmp2, true);
  cnt = 0;
  for (int i = 0; i < T.dim; i++) {
    mult = vtmp1[i] / v1L[i];
    cnt += abs(mult - lmd1) > test_eps;
  }
  printf("%d\n", cnt);

  printf("λ2 = %.12f\n", lmd2);
  vcopy(vtmp1, v2R, T.dim);
  T.product(temperature, vtmp1, vtmp2);
  cnt = 0;
  for (int i = 0; i < T.dim; i++) {
    mult = vtmp1[i] / v2R[i];
    cnt += abs(mult - lmd2) > test_eps && v2R[i] > test_eps;
  }
  printf("%d\n", cnt);

  printf("v1R・v1L = %.15f\n", dot(v1R, v1L, T.dim));
  printf("v2R・v1L = %.15f\n", dot(v2R, v1L, T.dim));


  for (int i = 0; i < T.dim; i++)
    printf("%.12f \n", v2R[i]);

  END();
}