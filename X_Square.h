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

class X_Square {
 public:
  double* Js;
  int M, dim, dim4, J_num;
  double EPS;

  X_Square(double* Js, int M, double EPS) : Js(Js), M(M), EPS(EPS), J_num(4) {
    dim = 1 << M;
    dim4 = dim << 2;

    double sum = 0.0;
    for (int i = 0; i < J_num; i++)
      sum += Js[i];
    for (int i = 0; i < J_num; i++)
      Js[i] /= sum * 2;
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
    double w3[2] = {1.0, exp(-2.0 * Js[3] / temperature)};

    // U1
    for (int s = 0; s < dim4; s++)
      vtmp[s] = 0.0;
    for (int s = 0; s < dim; s++) {
      int smm = s << 1;
      int spm = smm | 1, smp = smm | (dim << 1);
      int spp = spm | smp;

      int s_s[4] = {smm, smp, spm, spp};
      int sgm_ld[4] = {0, 0, 1, 1}, sgm_lu[4] = {0, 1, 0, 1};
      int sgm_rd = s & 1, sgm_rcd = (s >> 1) & 1, sgm_rcu = (s >> (M - 2)) & 1,
          sgm_ru = (s >> (M - 1)) & 1;

      for (int k = 0; k < 4; k++)
        vtmp[s_s[k]] += v[s] * w1[sgm_ld[k] ^ sgm_rd] * w2[sgm_ld[k] ^ sgm_ru] *
                        w3[sgm_ld[k] ^ sgm_rcd] * w1[sgm_lu[k] ^ sgm_ru] *
                        w2[sgm_lu[k] ^ sgm_rcu] * w3[sgm_lu[k] ^ sgm_rd];
    }
    for (int s = 0; s < dim4; s++)
      v[s] = vtmp[s];

    // U2 ~ U(M-1)
    for (int i = 1; i < M - 1; i++) {
      for (int s = 0; s < dim4; s++)
        vtmp[s] = 0.0;
      for (int s = 0; s < dim4; s++) {
        int sm = s & ~(1 << i), sp = s | (1 << i);

        int s_s[2] = {sm, sp};
        int sgm_l[2] = {0, 1};
        int sgm_rd = (s >> i) & 1, sgm_rc = (s >> (i + 1)) & 1,
            sgm_ru = (s >> (i + 2)) & 1;

        for (int k = 0; k < 2; k++)
          vtmp[s_s[k]] += v[s] * w1[sgm_l[k] ^ sgm_rc] * w2[sgm_l[k] ^ sgm_rd] *
                          w3[sgm_l[k] ^ sgm_ru];
      }
      for (int s = 0; s < dim4; s++)
        v[s] = vtmp[s];
    }

    // UM
    for (int s = 0; s < dim; s++)
      vtmp[s] = 0.0;
    for (int s = 0; s < dim; s++) {
      int mask_lu = dim >> 1;

      int smm = s | ((s & mask_lu) << 2);
      smm &= ~mask_lu;

      int smp = smm | dim, spm = smm | mask_lu;
      int spp = smp | spm;

      int s_s[4] = {smm, smp, spm, spp};
      for (int k = 0; k < 4; k++)
        vtmp[s] += v[s_s[k]];
    }
    for (int s = 0; s < dim; s++)
      v[s] = vtmp[s];
  }


  // Uv^T -> v
  void product_U_T(double temperature, double* v, double* vtmp) {
    double w1[2] = {1.0, exp(-2.0 * Js[1] / temperature)};
    double w2[2] = {1.0, exp(-2.0 * Js[2] / temperature)};
    double w3[2] = {1.0, exp(-2.0 * Js[3] / temperature)};

    for (int s = 0; s < dim4; s++)
      vtmp[s] = 0.0;
    for (int s = 0; s < dim; s++) {
      int mask_lu = dim >> 1;

      int smm = s | ((s & mask_lu) << 2);
      smm &= ~mask_lu;

      int smp = smm | dim, spm = smm | mask_lu;
      int spp = smp | spm;

      int s_s[4] = {smm, smp, spm, spp};
      for (int k = 0; k < 4; k++)
        vtmp[s_s[k]] += v[s];
    }
    for (int s = 0; s < dim4; s++)
      v[s] = vtmp[s];


    for (int i = M - 2; i >= 1; i--) {
      for (int s = 0; s < dim4; s++)
        vtmp[s] = 0.0;
      for (int s = 0; s < dim4; s++) {
        int sm = s & ~(1 << i), sp = s | (1 << i);

        int s_s[2] = {sm, sp};
        int sgm_l[2] = {0, 1};
        int sgm_rd = (s >> i) & 1, sgm_rc = (s >> (i + 1)) & 1,
            sgm_ru = (s >> (i + 2)) & 1;

        for (int k = 0; k < 2; k++)
          vtmp[s] += v[s_s[k]] * w1[sgm_l[k] ^ sgm_rc] * w2[sgm_l[k] ^ sgm_rd] *
                     w3[sgm_l[k] ^ sgm_ru];
      }
      for (int s = 0; s < dim4; s++)
        v[s] = vtmp[s];
    }


    for (int s = 0; s < dim; s++)
      vtmp[s] = 0.0;
    for (int s = 0; s < dim; s++) {
      int smm = s << 1;
      int spm = smm | 1, smp = smm | (dim << 1);
      int spp = spm | smp;

      int s_s[4] = {smm, smp, spm, spp};
      int sgm_ld[4] = {0, 0, 1, 1}, sgm_lu[4] = {0, 1, 0, 1};
      int sgm_rd = s & 1, sgm_rcd = (s >> 1) & 1, sgm_rcu = (s >> (M - 2)) & 1,
          sgm_ru = (s >> (M - 1)) & 1;

      for (int k = 0; k < 4; k++)
        vtmp[s] += v[s_s[k]] * w1[sgm_ld[k] ^ sgm_rd] * w2[sgm_ld[k] ^ sgm_ru] *
                   w3[sgm_ld[k] ^ sgm_rcd] * w1[sgm_lu[k] ^ sgm_ru] *
                   w2[sgm_lu[k] ^ sgm_rcu] * w3[sgm_lu[k] ^ sgm_rd];
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
