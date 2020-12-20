#include <vector>
#include "Xsquare.hpp"

int main() {
  START();
  vector<vector<double>> Js_s = {{0,0,1,500}};
  int M_start = 3, M_end = 15;

  FILE* fp;
  for (auto Js_vec : Js_s) {
    chrono_start = chrono::system_clock::now();
    for (int M = M_start; M <= M_end; M += 2) {
      int M1 = M, M2 = M1 + 2;
      double Js[4] = {Js_vec[0], Js_vec[1], Js_vec[2], Js_vec[3]};
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

        fp = fopen("./out/Xsquare/Js.txt", "a");
        for (int i = 0; i < 4; i++)
          fprintf(fp, "%.12f%c", X1.Js[i], (i == 3 ? '\n' : ' '));
        fclose(fp);

        fp = fopen("./out/Xsquare/Js_raw.txt", "a");
        for (int i = 0; i < 4; i++)
          fprintf(fp, "%.2f%c", Js[i], (i == 3 ? '\n' : ' '));
        fclose(fp);
      }

      fp = fopen("./out/Xsquare/Tc.txt", "a");
      fprintf(fp, "%.12f%c", Tc, (M == M_end ? '\n' : ' '));
      fclose(fp);
    }
    END();
  }
}