#include <chrono>
#include "Triangular.h"

int main() {
  chrono::system_clock::time_point start, end;
  start = chrono::system_clock::now();

  int M = 14;
  double Js[3] = {5.0, 1.0, 0.0};
  double temperature = 0.4;
  double EPS = 1e-12;

  Triangular T(Js, M, EPS);
  double* vo = alloc_dvector(T.dim);
  double* vn = alloc_dvector(T.dim);
  double* vtmp1 = alloc_dvector(T.dim2);
  double* vtmp2 = alloc_dvector(T.dim2);
  double* v1R = alloc_dvector(T.dim);
  double* v1L = alloc_dvector(T.dim);


  double lmd1 = T.power1_R(temperature, vo, vn, v1R, vtmp1, vtmp2);
  T.power1_L(temperature, vo, vn, v1L, vtmp1, vtmp2);
  double lmd2 = T.power2(temperature, lmd1, vo, vn, v1R, v1L, vtmp1, vtmp2);
  double xi = T.calc_xi(temperature, vo, vn, vtmp1, vtmp2, v1R, v1L);

  printf("λ1 = %.12f\n", lmd1);
  printf("λ2 = %.12f\n", lmd2);
  printf("ξ = %.12f\n", xi);

  end = chrono::system_clock::now();
  double time = static_cast<double>(
      chrono::duration_cast<chrono::microseconds>(end - start).count() /
      1000000.0);
  printf("time %lf[s]\n", time);
}