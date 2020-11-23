#include <chrono>
#include "Triangular.h"

int main() {
  chrono::system_clock::time_point start, end;
  start = chrono::system_clock::now();

  int M1 = 12;
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

  printf("%.4f\n", T1.exact_Tc());

  Triangular T2(Js, M2, EPS);
  double* vo_2 = alloc_dvector(T2.dim);
  double* vn_2 = alloc_dvector(T2.dim);
  double* vtmp1_2 = alloc_dvector(T2.dim2);
  double* vtmp2_2 = alloc_dvector(T2.dim2);
  double* v1R_2 = alloc_dvector(T2.dim);
  double* v1L_2 = alloc_dvector(T2.dim);

  int N = 20;
  double Tmin = 0.5, Tmax = 0.65;
  double dT = (Tmax - Tmin) / N;
  double* x = alloc_dvector(N + 1);
  double* y1 = alloc_dvector(N + 1);
  double* y2 = alloc_dvector(N + 1);
  for (int i = 0; i <= N; i++) {
    x[i] = Tmin + dT * i;
    double xi_1 = T1.calc_xi(x[i], vo_1, vn_1, vtmp1_1, vtmp2_1, v1R_1, v1L_1);
    double xi_2 = T2.calc_xi(x[i], vo_2, vn_2, vtmp1_2, vtmp2_2, v1R_2, v1L_2);
    y1[i] = M1 / xi_1;
    y2[i] = M2 / xi_2;
  }

  end = chrono::system_clock::now();
  double time = static_cast<double>(
      chrono::duration_cast<chrono::microseconds>(end - start).count() /
      1000000.0);
  printf("time %lf[s]\n", time);

  /* gnuplot */
  FILE* gp;
  gp = popen("gnuplot -persist", "w");
  fprintf(gp, "set xlabel \"T\"\n");
  fprintf(gp, "set ylabel \"M/xi\"\n");
  fprintf(gp, "set key left top\n");
  fprintf(gp,
          "plot "
          "'-' title \"M = %d\", "
          "'-' title \"M = %d\"\n",
          M1, M2);
  for (int i = 0; i <= N; i++) {
    fprintf(gp, "%.10f\t%.10f\n", x[i], y1[i]);
  }
  fprintf(gp, "e\n");
  for (int i = 0; i <= N; i++) {
    fprintf(gp, "%.10f\t%.10f\n", x[i], y2[i]);
  }
  fprintf(gp, "e\n");
  pclose(gp);
}