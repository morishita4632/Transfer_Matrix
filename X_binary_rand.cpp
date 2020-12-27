#include <vector>
#include "Xsquare.hpp"

int main() {
  int seed = now();
  printf("seed = %ld\n", seed);
  START(seed);

  int num = 10;
  int Js_range[4] = {1,100,1,1};
  
  int M_start = 3, M_end = 15;
  char buff_Js[128], buff_Js_raw[128], buff_Tc[128];
  sprintf (buff_Js, "./out/Xsquare/%d_%d_%d_%d_Js.txt", Js_range[0], Js_range[1], Js_range[2], Js_range[3]);
  sprintf (buff_Js_raw, "./out/Xsquare/%d_%d_%d_%d_Js_raw.txt", Js_range[0], Js_range[1], Js_range[2], Js_range[3]);
  sprintf (buff_Tc, "./out/Xsquare/%d_%d_%d_%d_Tc.txt", Js_range[0], Js_range[1], Js_range[2], Js_range[3]);

  FILE* fp;
  while (num-- > 0) {
    chrono_start = chrono::system_clock::now();
    double Js[4] = {rand01()*Js_range[0], rand01()*Js_range[1], rand01()*Js_range[2], rand01()*Js_range[3]};
    for (int M = M_start; M <= M_end; M += 2) {
      int M1 = M, M2 = M1 + 2;
      double EPS = 1e-12;
      double L = 0.1, R = 0.8;

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
        ((M1 / xi1 > M2 / xi2) ? l : r) = c;
      }

      double Tc = (l + r) / 2.0;

      if (M == M_start) {
        for (int i = 0; i < 4; i++)
          printf("%.2f%c", Js[i], (i == 3 ? '\n' : ' '));

        fp = fopen(buff_Js, "a");
        for (int i = 0; i < 4; i++)
          fprintf(fp, "%.12f%c", X1.Js[i], (i == 3 ? '\n' : ' '));
        fclose(fp);
        fp = fopen(buff_Js_raw, "a");
        for (int i = 0; i < 4; i++)
          fprintf(fp, "%.2f%c", Js[i], (i == 3 ? '\n' : ' '));
        fclose(fp);
      }
      fp = fopen(buff_Tc, "a");
      fprintf(fp, "%.12f%c", Tc, (M == M_end ? '\n' : ' '));
      fclose(fp);
    }
    END();
  }
}