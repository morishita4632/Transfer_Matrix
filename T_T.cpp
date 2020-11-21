#include "Triangular.h"

// やるときはTriangular.hのコンストラクタ内の規格化をオフにする

int main() {
  int M = 2;
  double Js[3] = {1.0, 1.0, 1.0};
  double temperature = 1.0;
  double EPS = 1e-12;

  Triangular T(Js, M, EPS);
  double* v = alloc_dvector(T.dim2);
  double* vtmp = alloc_dvector(T.dim2);

  for (int i = 0; i < T.dim; i++) {
    for (int s = 0; s < T.dim; s++)
      v[s] = (i == s);
    T.product(temperature, v, vtmp);

    for (int s = 0; s < T.dim; s++) {
      int val = log(v[s]) + 2;
      printf("%+3d%c", val, (s == T.dim - 1 ? '\n' : ' '));
    }
  }
}