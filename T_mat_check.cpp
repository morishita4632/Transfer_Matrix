#include "Triangular.h"

// TM自体が正しいか確かめる

void print_mat(double** mat, int size){
  for(int i=0; i<size; i++){
    for(int j=0; j<size; j++){
      printf("%10.9f%c", mat[i][j], (j==size-1)?'\n':' ');
    }
  }
}

int main() {
  START(0);

  int M = 4;
  double Js[3] = {rand01(), rand01(), rand01()};
  double EPS = 1e-12;
  double temperature = 0.6;

  Triangular T(Js, M, EPS);
  double** mat_exact = alloc_dmatrix(T.dim, T.dim);
  double** mat_calc = alloc_dmatrix(T.dim, T.dim);
  double* vtmp1 = alloc_dvector(T.dimt);
  double* vtmp2 = alloc_dvector(T.dimt);
  T.fill_T(temperature, mat_exact, vtmp1, vtmp2);
  T.fill_T_bruteforce(temperature, mat_calc);

  print_mat(mat_exact, T.dim);
  putchar('\n');
  print_mat(mat_calc, T.dim);

  END();
}