#include "Triangular.hpp"

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
  double Js[3] = {1, 1, 1};
  double EPS = 1e-12;
  double temperature = 1;

  Triangular T(Js, M, EPS);
  double** mat_exact = alloc_dmatrix(T.dim, T.dim);
  double* vtmp1 = alloc_dvector(T.dimt);
  double* vtmp2 = alloc_dvector(T.dimt);
  T.fill_T_bruteforce(temperature, mat_exact);

  for (int i = 0; i < T.dim; i++) {
    for (int j = 0; j < T.dim; j++) {
      mat_exact[i][j] = log(mat_exact[i][j]);
      mat_exact[i][j] /= 1.0 / 6.0;
      mat_exact[i][j] /= -2;
    }
  }

  print_mat(mat_exact, T.dim);

  END();
}