#include <math.h>
#include <chrono>
#include <random>
#include <iostream>
using namespace std;

chrono::system_clock::time_point chrono_start, chrono_end;
mt19937 mt(0);

static inline double* alloc_dvector(int n) {
  double* vec;
  vec = (double*)calloc(n, sizeof(double));
  if (vec == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_dvector\n");
    exit(1);
  }
  return vec;
}

/* allocate m x n column-major matrix of double */
static inline double **alloc_dmatrix(int m, int n) {
  int i;
  double **mat;
  mat = (double**)malloc((size_t)(n * sizeof(double*)));
  if (mat == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_dmatrix\n");
    exit(1);
  }
  mat[0] = (double*)calloc(m * n, sizeof(double));
  if (mat[0] == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_dmatrix\n");
    exit(1);
  }
  for (i = 1; i < n; ++i) mat[i] = mat[i-1] + m;
  return mat;
}

static inline double rand01() {
  return (double)(mt()) / (double)(mt.max());
}

static inline double dot(const double* v1, const double* v2, int dim) {
  double res = 0.0;
  for (int i = 0; i < dim; i++) {
    res += v1[i] * v2[i];
  }
  return res;
}

static inline void vcopy(double* v1, const double* v2, int dim) {
  for (int i = 0; i < dim; i++) {
    v1[i] = v2[i];
  }
}

static inline void vzero(double* v, int dim) {
  for (int i = 0; i < dim; i++) {
    v[i] = 0.0;
  }
}

// Measure time & Random initialize
static inline void START(int seed = 0) {
  mt.seed(seed);
  chrono_start = chrono::system_clock::now();
}

static inline void END() {
  chrono_end = chrono::system_clock::now();
  double time = static_cast<double>(
      chrono::duration_cast<chrono::microseconds>(chrono_end - chrono_start)
          .count() /
      1000000.0);
  printf("time %.2lf[s],  %.2lf[m]\n", time, time / 60.0);
}