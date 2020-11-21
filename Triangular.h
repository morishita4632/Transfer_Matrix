#include <math.h>
#include <iostream>
using namespace std;

/* allocate vector of double */
static inline double* alloc_dvector(int n) {
  double* vec;
  vec = (double*)calloc(n, sizeof(double));
  if (vec == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_dvector\n");
    exit(1);
  }
  return vec;
}

double rand01() {
  return (double)rand() / (double)RAND_MAX;
}

class Triangular {
 public:
  double* Js;
  int M, dim, dim2, J_num;
  double EPS;

  Triangular(double* Js, int M, double EPS) : Js(Js), M(M), EPS(EPS), J_num(3) {
    dim = 1 << M;
    dim2 = dim << 1;

    double sum = 0.0;
    for (int i = 0; i < J_num; i++)
      sum += Js[i];
    for (int i = 0; i < J_num; i++)
      Js[i] /= sum * 2;
  }

  double exact_Tc() {
    double l = 0.10, r = 0.61, c;
    while (r - l > EPS) {
      c = (l + r) / 2.0;

      double f = 0.0;
      for (int i = 0; i < J_num; i++) {
        int j = (i + 1) % J_num;
        f += exp(-2.0 * (Js[i] + Js[j]) / c);
      }

      (f > 1.0 ? r : l) = c;
    }
    return (l + r) / 2.0;
  }

  // Dv -> v
  void product_D(double temperature, double* v) {
    double w0[2] = {1.0, exp(-2.0 * Js[0] / temperature)};
    for (int s = 0; s < dim; s++) {
      for (int i = 0; i < M; i++) {
        int j = (i + 1) % M;
        int sigma_i = (s >> i) & 1;
        int sigma_j = (s >> j) & 1;
        v[s] = v[s] * w0[sigma_i ^ sigma_j];
      }
    }
  }

  // Uv -> v
  void product_U(double temperature, double* v, double* vtmp) {
    double w1[2] = {1.0, exp(-2.0 * Js[1] / temperature)};
    double w2[2] = {1.0, exp(-2.0 * Js[2] / temperature)};

    // U1
    for (int s = 0; s < dim2; s++)
      vtmp[s] = 0.0;
    for (int s = 0; s < dim; s++) {
      int sm = s << 1, sp = (s << 1) | 1;
      int sigma_d = s & 1, sigma_u = (s >> (M - 1)) & 1;
      vtmp[sp] += v[s] * w1[1 ^ sigma_d] * w2[1 ^ sigma_u];
      vtmp[sm] += v[s] * w1[0 ^ sigma_d] * w2[0 ^ sigma_u];
    }
    for (int s = 0; s < dim2; s++)
      v[s] = vtmp[s];

    // U2 ~ UM
    for (int i = 1; i < M; i++) {
      for (int s = 0; s < dim2; s++)
        vtmp[s] = 0.0;
      for (int s = 0; s < dim2; s++) {
        int sm = s & ~(1 << i), sp = s | (1 << i);
        int sigma_d = (s >> i) & 1, sigma_u = (s >> (i + 1)) & 1;
        vtmp[sp] += v[s] * w1[1 ^ sigma_u] * w2[1 ^ sigma_d];
        vtmp[sm] += v[s] * w1[0 ^ sigma_u] * w2[0 ^ sigma_d];
      }
      for (int s = 0; s < dim2; s++)
        v[s] = vtmp[s];
    }

    // U(M+1)
    for (int s = 0; s < dim; s++)
      vtmp[s] = 0.0;
    for (int s = 0; s < dim; s++) {
      int sm = s, sp = s | dim;
      vtmp[s] += v[sm];
      vtmp[s] += v[sp];
    }
    for (int s = 0; s < dim; s++)
      v[s] = vtmp[s];
  }

  // U^T v -> v
  void product_U_T(double temperature, double* v, double* vtmp) {
    double w1[2] = {1.0, exp(-2.0 * Js[1] / temperature)};
    double w2[2] = {1.0, exp(-2.0 * Js[2] / temperature)};

    for (int s = 0; s < dim2; s++)
      vtmp[s] = 0.0;
    for (int s = 0; s < dim; s++) {
      int sm = s, sp = s | dim;
      vtmp[sm] += v[s];
      vtmp[sp] += v[s];
    }
    for (int s = 0; s < dim2; s++)
      v[s] = vtmp[s];

    for (int i = M - 1; i >= 1; i--) {
      for (int s = 0; s < dim2; s++)
        vtmp[s] = 0.0;
      for (int s = 0; s < dim2; s++) {
        int sm = s & ~(1 << i), sp = s | (1 << i);
        int sigma_d = (s >> i) & 1, sigma_u = (s >> (i + 1)) & 1;
        vtmp[s] += v[sp] * w1[1 ^ sigma_u] * w2[1 ^ sigma_d];
        vtmp[s] += v[sm] * w1[0 ^ sigma_u] * w2[0 ^ sigma_d];
      }
      for (int s = 0; s < dim2; s++)
        v[s] = vtmp[s];
    }

    for (int s = 0; s < dim; s++)
      vtmp[s] = 0.0;
    for (int s = 0; s < dim; s++) {
      int sm = s << 1, sp = (s << 1) | 1;
      int sigma_d = s & 1, sigma_u = (s >> (M - 1)) & 1;
      vtmp[s] += v[sp] * w1[1 ^ sigma_d] * w2[1 ^ sigma_u];
      vtmp[s] += v[sm] * w1[0 ^ sigma_d] * w2[0 ^ sigma_u];
    }
    for (int s = 0; s < dim; s++)
      v[s] = vtmp[s];
  }

  // Tv -> v
  void product(double temperature, double* v, double* vtmp,
               bool transpose = false) {
    if (transpose) {
      product_D(temperature, v);
      product_U_T(temperature, v, vtmp);
    } else {
      product_U(temperature, v, vtmp);
      product_D(temperature, v);
    }
  }

  // return λ1, fill v1R
  double power1_R(double temperature, double* vo, double* vn, double* v1R,
                  double* vtmp1, double* vtmp2) {
    for (int s = 0; s < dim; s++)
      vo[s] = rand01();
    double lmd_o = 0.0, lmd_n = 0.0;
    do {
      lmd_o = lmd_n;
      for (int s = 0; s < dim; s++)
        vtmp1[s] = vo[s];
      product(temperature, vtmp1, vtmp2);
      for (int s = 0; s < dim; s++)
        vn[s] = vtmp1[s];

      double numer = 0.0, denom = 0.0;
      for (int s = 0; s < dim; s++) {
        numer += vn[s] * vn[s];
        denom += vn[s] * vo[s];
      }
      lmd_n = numer / denom;
      for (int s = 0; s < dim; s++)
        vo[s] = vn[s] / sqrt(numer);
    } while (abs(lmd_n - lmd_o) > EPS);

    for (int s = 0; s < dim; s++)
      v1R[s] = vo[s];
    return lmd_n;
  }

  // fill v1L
  void power1_L(double temperature, double* vo, double* vn, double* v1L,
                double* vtmp1, double* vtmp2) {
    for (int s = 0; s < dim; s++)
      vo[s] = rand01();
    double lmd_o = 0.0, lmd_n = 0.0;
    do {
      lmd_o = lmd_n;
      for (int s = 0; s < dim; s++)
        vtmp1[s] = vo[s];
      product(temperature, vtmp1, vtmp2, true);
      for (int s = 0; s < dim; s++)
        vn[s] = vtmp1[s];

      double numer = 0.0, denom = 0.0;
      for (int s = 0; s < dim; s++) {
        numer += vn[s] * vn[s];
        denom += vn[s] * vo[s];
      }
      lmd_n = numer / denom;
      for (int s = 0; s < dim; s++)
        vo[s] = vn[s] / sqrt(numer);
    } while (abs(lmd_n - lmd_o) > EPS);

    for (int s = 0; s < dim; s++)
      v1L[s] = vo[s];
  }

  // λ・vR・vL -> vtmp
  void decrease_term(double lmd1, const double* v1R, const double* v1L,
                     const double* v, double* vtmp) {
    double v1L_inner_v = 0.0;
    for (int s = 0; s < dim; s++)
      v1L_inner_v += v1L[s] * v[s];
    for (int s = 0; s < dim; s++)
      vtmp[s] = lmd1 * v1L_inner_v * v1R[s];
  }

  // return λ2
  double power2(double temperature, double lmd1, double* vo, double* vn,
                const double* v1R, const double* v1L, double* vtmp1,
                double* vtmp2) {
    for (int s = 0; s < dim; s++)
      vo[s] = rand01();
    double lmd_o = 0.0, lmd_n = 0.0;
    do {
      lmd_o = lmd_n;
      for (int s = 0; s < dim; s++)
        vtmp1[s] = vo[s];
      product(temperature, vtmp1, vtmp2);
      for (int s = 0; s < dim; s++)
        vn[s] = vtmp1[s];
      decrease_term(lmd1, v1R, v1L, vo, vtmp1);
      for (int s = 0; s < dim; s++)
        vn[s] -= vtmp1[s];

      double numer = 0.0, denom = 0.0;
      for (int s = 0; s < dim; s++) {
        numer += vn[s] * vn[s];
        denom += vn[s] * vo[s];
      }
      lmd_n = numer / denom;
      for (int s = 0; s < dim; s++)
        vo[s] = vn[s] / sqrt(numer);
    } while (abs(lmd_n - lmd_o) > EPS);

    return lmd_n;
  }

  double calc_xi(double temperature, double* vo, double* vn, double* vtmp1,
                 double* vtmp2, double* v1R, double* v1L) {
    double lmd1 = power1_R(temperature, vo, vn, v1R, vtmp1, vtmp2);
    power1_L(temperature, vo, vn, v1L, vtmp1, vtmp2);
    double lmd2 = power2(temperature, lmd1, vo, vn, v1R, v1L, vtmp1, vtmp2);
    return 1.0 / abs(log(lmd2 / lmd1));
  }
};
