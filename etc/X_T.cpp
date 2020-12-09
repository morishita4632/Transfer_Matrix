#include "X_Square.hpp"

int main() {
  int M = 3;
  double Js[4] = {1.0, 1.0, 1.0, 0.0};
  double temperature = 1.0;
  double EPS = 1e-12;

  X_Square X(Js, M, EPS);
  double* v = alloc_dvector(X.dim4);
  double* vtmp = alloc_dvector(X.dim4);

  for (int i = 0; i < X.dim; i++) {
    for (int s = 0; s < X.dim; s++)
      v[s] = (i == s);
    X.product(temperature, v, vtmp);

    for (int s = 0; s < X.dim; s++) {
      // int val = log(v[s]) + 2;
      // if (val == -2147483648)
      //   printf("   %c", (s == X.dim - 1 ? '\n' : ' '));
      // else
      //   printf("%+3d%c", val, (s == X.dim - 1 ? '\n' : ' '));
      printf("%f%c", v[s], (s == X.dim - 1 ? '\n' : ' '));
    }
  }
}
