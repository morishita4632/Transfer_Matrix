#include "Triangular.hpp"

// TMを愚直解と比較する

void print_mat(double** mat, int size) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      printf("%10.6f%c", mat[i][j], (j == size - 1) ? '\n' : ' ');
    }
  }
}

int main() {
  START(1);

  double TEST_EPS = 1e-8;

  int M = 6;
  double temperature = 0.6;
  double EPS = 1e-12;


  int NG = 0;
  for (int n = 0; n < 1000; n++) {
    double Js[3] = {rand01(), rand01(), rand01()};

    Triangular T(Js, M, EPS);
    double** mat_exact = alloc_dmatrix(T.dim, T.dim);
    double** mat_calc = alloc_dmatrix(T.dim, T.dim);
    double* vtmp1 = alloc_dvector(T.dimt);
    double* vtmp2 = alloc_dvector(T.dimt);
    T.fill_T_bruteforce(temperature, mat_exact);
    T.fill_T(temperature, mat_calc, vtmp1, vtmp2);

    for (int i = 0; i < T.dim; i++) {
      for (int j = 0; j < T.dim; j++) {
        NG += abs(mat_calc[i][j] - mat_exact[i][j]) > TEST_EPS;
      }
    }
  }
  cout << NG << endl;

  END();
}