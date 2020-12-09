#include "Triangular.hpp"

int main() {
  START();

  int M1 = 9;
  int M2 = M1 + 1;
  double Js[3] = {7.0, 1.0, 1.0};
  double EPS = 1e-12;

  Triangular T1(Js, M1, EPS);
  double* vo_1 = alloc_dvector(T1.dim);
  double* vn_1 = alloc_dvector(T1.dim);
  double* vtmp1_1 = alloc_dvector(T1.dimt);
  double* vtmp2_1 = alloc_dvector(T1.dimt);
  double* v1R_1 = alloc_dvector(T1.dim);
  double* v1L_1 = alloc_dvector(T1.dim);
  double* v2R_1 = alloc_dvector(T1.dim);

  double Tc_ex = T1.exact_Tc();
  printf("%.4f\n", Tc_ex);

  Triangular T2(Js, M2, EPS);
  double* vo_2 = alloc_dvector(T2.dim);
  double* vn_2 = alloc_dvector(T2.dim);
  double* vtmp1_2 = alloc_dvector(T2.dimt);
  double* vtmp2_2 = alloc_dvector(T2.dimt);
  double* v1R_2 = alloc_dvector(T2.dim);
  double* v1L_2 = alloc_dvector(T2.dim);
  double* v2R_2 = alloc_dvector(T2.dim);

  int N = 20;
  double Tmin = 0.4, Tmax = 0.6;
  double dT = (Tmax - Tmin) / N;
  double* x = alloc_dvector(N + 1);
  double* y1 = alloc_dvector(N + 1);
  double* y2 = alloc_dvector(N + 1);
  for (int i = 0; i <= N; i++) {
    x[i] = Tmin + dT * i;
    double xi_1 =
        T1.calc_xi(x[i], vo_1, vn_1, vtmp1_1, vtmp2_1, v1R_1, v1L_1, v2R_1);
    double xi_2 =
        T2.calc_xi(x[i], vo_2, vn_2, vtmp1_2, vtmp2_2, v1R_2, v1L_2, v2R_2);
    y1[i] = M1 / xi_1;
    y2[i] = M2 / xi_2;
  }

  END();

  double ymin = 0.0, ymax = 4.0;
  /* gnuplot */
  FILE* gp;
  gp = popen("gnuplot -persist", "w");
  fprintf(gp, "set xlabel \"T\"\n");
  fprintf(gp, "set ylabel \"M/xi\"\n");
  fprintf(gp, "set yrange [%f:%f]\n", ymin, ymax);
  fprintf(gp, "set key left top\n");
  fprintf(gp, "set parametric\n");
  fprintf(gp, "set trange [%f:%f]\n", ymin, ymax);
  fprintf(gp,
          "plot "
          "'-' title \"M = %d\", "
          "'-' title \"M = %d\", "
          "%f, t with lines notitle\n",
          M1, M2, Tc_ex);
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