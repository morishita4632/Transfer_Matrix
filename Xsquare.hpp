#include "utility.hpp"

class Xsquare {
 public:
  double* Js;
  int M, dim, dimt, J_num;
  double EPS;

  Xsquare(double* Js, int M, double EPS) : M(M), EPS(EPS), J_num(4) {
    dim = 1 << M;
    dimt = dim << 2;

    this->Js = alloc_dvector(J_num);
    vcopy(this->Js, Js, J_num);
    double sum = 0.0;
    for (int i = 0; i < J_num; i++)
      sum += this->Js[i];
    for (int i = 0; i < J_num; i++)
      this->Js[i] /= sum * 2;
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
    vzero(vtmp, dimt);
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
    vcopy(v, vtmp, dimt);

    // U2 ~ U(M-1)
    for (int i = 1; i < M - 1; i++) {
      vzero(vtmp, dimt);
      for (int s = 0; s < dimt; s++) {
        int sm = s & ~(1 << i), sp = s | (1 << i);

        int s_s[2] = {sm, sp};
        int sgm_l[2] = {0, 1};
        int sgm_rd = (s >> i) & 1, sgm_rc = (s >> (i + 1)) & 1,
            sgm_ru = (s >> (i + 2)) & 1;

        for (int k = 0; k < 2; k++)
          vtmp[s_s[k]] += v[s] * w1[sgm_l[k] ^ sgm_rc] * w2[sgm_l[k] ^ sgm_rd] *
                          w3[sgm_l[k] ^ sgm_ru];
      }
      vcopy(v, vtmp, dimt);
    }

    // UM
    vzero(vtmp, dim);
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
    vcopy(v, vtmp, dim);
  }


  // Uv^T -> v
  void product_U_T(double temperature, double* v, double* vtmp) {
    double w1[2] = {1.0, exp(-2.0 * Js[1] / temperature)};
    double w2[2] = {1.0, exp(-2.0 * Js[2] / temperature)};
    double w3[2] = {1.0, exp(-2.0 * Js[3] / temperature)};

    vzero(vtmp, dimt);
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
    vcopy(v, vtmp, dimt);


    for (int i = M - 2; i >= 1; i--) {
      vzero(vtmp, dimt);
      for (int s = 0; s < dimt; s++) {
        int sm = s & ~(1 << i), sp = s | (1 << i);

        int s_s[2] = {sm, sp};
        int sgm_l[2] = {0, 1};
        int sgm_rd = (s >> i) & 1, sgm_rc = (s >> (i + 1)) & 1,
            sgm_ru = (s >> (i + 2)) & 1;

        for (int k = 0; k < 2; k++)
          vtmp[s] += v[s_s[k]] * w1[sgm_l[k] ^ sgm_rc] * w2[sgm_l[k] ^ sgm_rd] *
                     w3[sgm_l[k] ^ sgm_ru];
      }
      vcopy(v, vtmp, dimt);
    }

    vzero(vtmp, dim);
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
    vcopy(v, vtmp, dim);
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
      vcopy(vtmp1, vo, dim);
      product(temperature, vtmp1, vtmp2);
      vcopy(vn, vtmp1, dim);

      double numer = dot(vn, vn, dim), denom = dot(vn, vo, dim);
      lmd_n = numer / denom;

      for (int s = 0; s < dim; s++)
        vo[s] = vn[s] / sqrt(numer);
    } while (abs(lmd_n - lmd_o) > EPS);

    vcopy(v1R, vo, dim);
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
      vcopy(vtmp1, vo, dim);
      product(temperature, vtmp1, vtmp2, true);
      vcopy(vn, vtmp1, dim);

      double numer = dot(vn, vn, dim), denom = dot(vn, vo, dim);
      lmd_n = numer / denom;

      for (int s = 0; s < dim; s++)
        vo[s] = vn[s] / sqrt(numer);
    } while (abs(lmd_n - lmd_o) > EPS);

    vcopy(v1L, vo, dim);
  }

  // λ・vR・vL -> vtmp
  void decrease_term(double lmd1, const double* v1R, const double* v1L,
                     const double* v, double* vtmp) {
    double v1L_inner_v = dot(v1L, v, dim);
    for (int s = 0; s < dim; s++)
      vtmp[s] = lmd1 * v1L_inner_v * v1R[s];
  }

  // return λ2, fill v2R
  double power2(double temperature, double lmd1, double* vo, double* vn,
                const double* v1R, const double* v1L, double* v2R,
                double* vtmp1, double* vtmp2) {
    for (int s = 0; s < dim; s++)
      vo[s] = rand01();
    double lmd_o = 0.0, lmd_n = 0.0;
    do {
      lmd_o = lmd_n;
      vcopy(vtmp1, vo, dim);
      product(temperature, vtmp1, vtmp2);
      vcopy(vn, vtmp1, dim);
      decrease_term(lmd1, v1R, v1L, vo, vtmp1);
      for (int s = 0; s < dim; s++)
        vn[s] -= vtmp1[s];

      double numer = dot(vn, vn, dim), denom = dot(vn, vo, dim);
      lmd_n = numer / denom;

      for (int s = 0; s < dim; s++)
        vo[s] = vn[s] / sqrt(numer);
    } while (abs(lmd_n - lmd_o) > EPS);

    vcopy(v2R, vo, dim);
    return lmd_n;
  }

  // Return xi. If (v2R, v1L)!=0, return -1.0;
  double calc_xi(double temperature, double* vo, double* vn, double* vtmp1,
                 double* vtmp2, double* v1R, double* v1L, double* v2R) {
    double lmd1 = power1_R(temperature, vo, vn, v1R, vtmp1, vtmp2);
    power1_L(temperature, vo, vn, v1L, vtmp1, vtmp2);
    double inner = dot(v1R, v1L, dim);
    for (int i = 0; i < dim; i++)
      v1L[i] /= inner;
    double lmd2 =
        power2(temperature, lmd1, vo, vn, v1R, v1L, v2R, vtmp1, vtmp2);
    double xi = 1.0 / log(lmd1 / lmd2);
    return xi;
  }

  // For confirmation
  void fill_T(double temperature, double** mat, double* vtmp1, double* vtmp2) {
    for (int j = 0; j < dim; j++) {
      vzero(vtmp1, dimt);
      vtmp1[j] = 1.0;
      product(temperature, vtmp1, vtmp2);
      for (int i = 0; i < dim; i++)
        mat[i][j] = vtmp1[i];
    }
  }

  void fill_T_bruteforce(double temperature, double** mat) {
    double w[4][2] = {{1.0, exp(-2.0 * Js[0] / temperature)},
                      {1.0, exp(-2.0 * Js[1] / temperature)},
                      {1.0, exp(-2.0 * Js[2] / temperature)},
                      {1.0, exp(-2.0 * Js[3] / temperature)}};
    for (int s1 = 0; s1 < dim; s1++) {
      for (int s2 = 0; s2 < dim; s2++) {
        double temp = 1.0;
        for (int i = 0; i < M; i++) {
          int j = (i + M - 1) % M;
          int jj = (i + 1) % M;
          int sgm = (s1 >> i) & 1;
          int sgm_s[4] = {(s1 >> j) & 1, (s2 >> i) & 1, (s2 >> j) & 1,
                          (s2 >> jj) & 1};

          for (int k = 0; k < 4; k++) {
            temp *= w[k][sgm ^ sgm_s[k]];
          }
        }
        mat[s1][s2] = temp;
      }
    }
  }
};
