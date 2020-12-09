#include "Xsquare.hpp"

// TMの愚直解を一応確かめる
// anti-ferro の個数で表示。PBCに注意

void print_mat(double** mat, int size) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      printf("%d%c", (int)(mat[i][j] + 0.5), (j == size - 1) ? '\n' : ' ');
    }
  }
}

int main() {
  START(0);

  int M = 2;
  double Js[4] = {1, 1, 1, 1};
  double EPS = 1e-12;
  double temperature = 1;

  Xsquare X(Js, M, EPS);
  double** mat_exact = alloc_dmatrix(X.dim, X.dim);
  double* vtmp1 = alloc_dvector(X.dimt);
  double* vtmp2 = alloc_dvector(X.dimt);
  X.fill_T_bruteforce(temperature, mat_exact);

  for (int i = 0; i < X.dim; i++) {
    for (int j = 0; j < X.dim; j++) {
      mat_exact[i][j] = log(mat_exact[i][j]);
      mat_exact[i][j] /= 1.0 / 8.0;
      mat_exact[i][j] /= -2;
    }
  }

  print_mat(mat_exact, X.dim);
  END();
}