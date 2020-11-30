#include <math.h>
#include <chrono>
#include <iostream>
using namespace std;

static inline double* alloc_dvector(int n) {
  double* vec;
  vec = (double*)calloc(n, sizeof(double));
  if (vec == NULL) {
    fprintf(stderr, "Error: allocation failed in alloc_dvector\n");
    exit(1);
  }
  return vec;
}

static inline double rand01() {
  return (double)rand() / (double)RAND_MAX;
}

static inline double dot(const double* v1, const double* v2, int dim) {
  double res = 0.0;
#pragma omp parallel for
  for (int i = 0; i < dim; i++) {
    res += v1[i] * v2[i];
  }
  return res;
}

static inline void vcopy(double* v1, const double* v2, int dim) {
#pragma omp parallel for
  for (int i = 0; i < dim; i++) {
    v1[i] = v2[i];
  }
}

static inline void vzero(double* v, int dim) {
#pragma omp parallel for
  for (int i = 0; i < dim; i++) {
    v[i] = 0.0;
  }
}

// Measure time
chrono::system_clock::time_point chrono_start, chrono_end;

static inline void START() {
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