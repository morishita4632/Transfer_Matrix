#include "Triangular.h"

int main() {
  chrono::system_clock::time_point start, end;
  start = chrono::system_clock::now();

  int M = 6;
  double Js[3] = {7.0, 1.0, 1.0};
  double EPS = 1e-15;

  double temperature = 0.6;

  Triangular T(Js, M, EPS);
  double* vo = alloc_dvector(T.dim);
  double* vn = alloc_dvector(T.dim);
  double* vtmp1 = alloc_dvector(T.dim2);
  double* vtmp2 = alloc_dvector(T.dim2);
  double* v1R = alloc_dvector(T.dim);
  double* v1L = alloc_dvector(T.dim);
  double* v2R = alloc_dvector(T.dim);

  double lmd1 = T.power1_R(temperature, vo, vn, v1R, vtmp1, vtmp2);
  T.power1_L(temperature, vo, vn, v1L, vtmp1, vtmp2);
  double lmd2 =
      T.power2(temperature, lmd1, vo, vn, v1R, v1L, vtmp1, vtmp2, v2R);
  printf("%.9f, %.9f\n", lmd1, lmd2);

  end = chrono::system_clock::now();
  double time = static_cast<double>(
      chrono::duration_cast<chrono::microseconds>(end - start).count() /
      1000000.0);
  printf("time %lf[s]\n", time);
}