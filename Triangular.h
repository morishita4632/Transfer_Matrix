#include "utility.h"

class Triangular {
 public:
  double* Js;
  int M, dim, dimt, J_num;
  double EPS;

  Triangular(double* Js, int M, double EPS) : M(M), EPS(EPS), J_num(3) {
    dim = 1 << M;
    dimt = dim << 1;

    this->Js = alloc_dvector(J_num);
    vcopy(this->Js, Js, J_num);
    double sum = 0.0;
    for (int i = 0; i < J_num; i++)
      sum += this->Js[i];
    for (int i = 0; i < J_num; i++)
      this->Js[i] /= sum * 2;
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
    vzero(vtmp, dimt);
    for (int s = 0; s < dim; s++) {
      int sm = s << 1, sp = (s << 1) | 1;
      int sigma_d = s & 1, sigma_u = (s >> (M - 1)) & 1;
      vtmp[sp] += v[s] * w1[1 ^ sigma_d] * w2[1 ^ sigma_u];
      vtmp[sm] += v[s] * w1[0 ^ sigma_d] * w2[0 ^ sigma_u];
    }
    vcopy(v, vtmp, dimt);

    // U2 ~ UM
    for (int i = 1; i < M; i++) {
      vzero(vtmp, dimt);
      for (int s = 0; s < dimt; s++) {
        int sm = s & ~(1 << i), sp = s | (1 << i);
        int sigma_d = (s >> i) & 1, sigma_u = (s >> (i + 1)) & 1;
        vtmp[sp] += v[s] * w1[1 ^ sigma_u] * w2[1 ^ sigma_d];
        vtmp[sm] += v[s] * w1[0 ^ sigma_u] * w2[0 ^ sigma_d];
      }
      vcopy(v, vtmp, dimt);
    }

    // U(M+1)
    vzero(vtmp, dim);
    for (int s = 0; s < dim; s++) {
      int sm = s, sp = s | dim;
      vtmp[s] += v[sm];
      vtmp[s] += v[sp];
    }
    vcopy(v, vtmp, dim);
  }

  // Tv -> v
  void product(double temperature, double* v, double* vtmp) {
    product_U(temperature, v, vtmp);
    product_D(temperature, v);
  }

  // return λ1, fill v1
  double power1(double temperature, double* vo, double* vn, double* v1,
                double* vtmp1, double* vtmp2) {
    for (int s = 0; s < dim; s++)
      vo[s] = rand01();
    double lmd_o = 0.0, lmd_n = 0.0;
    do {
      lmd_o = lmd_n;
      vcopy(vtmp1, vo, dim);
      product(temperature, vtmp1, vtmp2);
      vcopy(vn, vtmp1, dim);

      double numer = dot(vn, vn, dim), denom = dot(vn, vo, dim);
      lmd_n = numer / denom;

      for (int s = 0; s < dim; s++)
        vo[s] = vn[s] / sqrt(numer);
    } while (abs(lmd_n - lmd_o) > EPS);

    vcopy(v1, vo, dim);
    return lmd_n;
  }

  // λ1・v1・v1T -> vtmp
  void decrease_term(double lmd1, const double* v1, const double* v,
                     double* vtmp) {
    double v1_inner_v = dot(v1, v, dim);
    for (int s = 0; s < dim; s++)
      vtmp[s] = lmd1 * v1_inner_v * v1[s];
  }

  // return λ2, fill v2
  double power2(double temperature, double lmd1, double* vo, double* vn,
                const double* v1, double* v2, double* vtmp1, double* vtmp2) {
    for (int s = 0; s < dim; s++)
      vo[s] = rand01();
    double lmd_o = 0.0, lmd_n = 0.0;
    do {
      lmd_o = lmd_n;
      vcopy(vtmp1, vo, dim);
      product(temperature, vtmp1, vtmp2);
      vcopy(vn, vtmp1, dim);
      decrease_term(lmd1, v1, vo, vtmp1);
      for (int s = 0; s < dim; s++)
        vn[s] -= vtmp1[s];

      double numer = dot(vn, vn, dim), denom = dot(vn, vo, dim);
      lmd_n = numer / denom;

      for (int s = 0; s < dim; s++)
        vo[s] = vn[s] / sqrt(numer);
    } while (abs(lmd_n - lmd_o) > EPS);

    vcopy(v2, vo, dim);
    return lmd_n;
  }

  // Return xi
  double calc_xi(double temperature, double* vo, double* vn, double* vtmp1,
                 double* vtmp2, double* v1, double* v2) {
    double lmd1 = power1(temperature, vo, vn, v1, vtmp1, vtmp2);
    double lmd2 = power2(temperature, lmd1, vo, vn, v1, v2, vtmp1, vtmp2);
    double xi = 1.0 / abs(log(lmd2 / lmd1));
    return xi;
  }
};
