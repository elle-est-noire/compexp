#include "cmatrix.h"
#include "mersenne_twister.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define I _Complex_I
double range = 100.0;

#define SIZE 5

int main(int argc, char **argv) {
  char *filename;
  FILE *fp;

  int m, n;
  double complex **a, **a2; /* matrix */
  double *w;        /* eigenvalues */
  double complex *work;     /* working area */

  int lwork, info;
  char jobz = 'V'; /* Compute eigenvalues and eigenvectors */
  char uplo = 'U'; /* Upper triangle of A is stored */

  init_genrand((unsigned long)time(NULL));
  n = SIZE;
  a = alloc_zmatrix(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < i; j++) {
      a[i][j] = (genrand_real3() * 2.0 * range - range) + (genrand_real3() * 2.0 * range - range) * I;
      a[j][i] = conj(a[i][j]);
    }
    a[i][i] = genrand_real3() * 2.0 * range - range;
  }

  printf("Matrix A:\n");
  fprint_zmatrix(stdout, n, n, a);

  a2 = alloc_zmatrix(n, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      a2[i][j] = a[i][j];

  /* perform eigenvalue decomposition */
  w = alloc_dvector(n);
  // lwork = 3 * n - 1;
  lwork = 6 * n;
  double rwork[3 * SIZE - 2];
  work = alloc_zvector(lwork);
  zheev_(&jobz, &uplo, &n, mat_ptr(a), &n, vec_ptr(w), vec_ptr(work), &lwork, rwork, &info);
  if (info != 0) {
    fprintf(stderr, "Error: LAPACK::dsyev failed\n");
    exit(1);
  }
  printf("Eigenvalues:\n");
  fprint_dvector(stdout, n, w);
  printf("Eigenvectors [each column represents each eigenvector]:\n");
  fprint_zmatrix(stdout, n, n, a);

  // double complex **U;
  // U = alloc_zmatrix(n, n);
  // for (int i = 0; i < n; i++) {
  //   for (int j = 0; j < n; j++) {
  //     U[i][j] = 0.0;
  //     for (int k = 0; k < n; k++) {
  //       U[i][j] += a[i][k] * a[j][k];
  //     }
  //   }
  // }
  // fprint_zmatrix(stdout, n, n, U);

  // double complex **V;
  // V = alloc_zmatrix(n, n);
  // for (int i = 0; i < n; i++) {
  //   for (int j = 0; j < n; j++) {
  //     U[i][j] = 0.0;
  //     for (int k = 0; k < n; k++) {
  //       U[i][j] += a[i][k] * a2[k][j];
  //     }
  //   }
  // }
  // for (int i = 0; i < n; i++) {
  //   for (int j = 0; j < n; j++) {
  //     V[i][j] = 0.0;
  //     for (int k = 0; k < n; k++) {
  //       V[i][j] += U[i][k] * U[j][k];
  //     }
  //   }
  // }
  // fprint_zmatrix(stdout, n, n, V);
  // for (int i = 0; i < n; i++) {
  //   w[i] *= w[i];
  // }
  // fprint_zvector(stdout, n, w);

  double complex **U;
  U = alloc_zmatrix(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      U[i][j] = 0.0;
      for (int k = 0; k < n; k++) {
        U[i][j] += a[i][k] * conj(a[j][k]);
      }
      if (j == i) U[i][j] -= 1.0;
    }
  }
  printf("Matrix (tU*)U-I:\n");
  fprint_zmatrix_e(stdout, n, n, U);

  double complex **V;
  V = alloc_zmatrix(n, n);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      U[i][j] = 0.0;
      for (int k = 0; k < n; k++) {
        U[i][j] += a[i][k] * a2[k][j];
      }
    }
  }
  printf("Matrix V:\n");
  fprint_zmatrix(stdout, n, n, U);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      V[i][j] = a[i][j] * w[i];
    }
  }
  printf("Matrix V':\n");
  fprint_zmatrix(stdout, n, n, V);

  double *diff;
  diff = alloc_dvector(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      diff[i] += cabs(U[i][j] - V[i][j]) * cabs(U[i][j] - V[i][j]);
    }
    diff[i] = sqrt(diff[i]);
  }
  printf("diffs:\n");
  fprint_dvector(stdout, n, diff);

  free_zmatrix(a);
  free_dvector(w);
  free_zvector(work);
  free_dvector(diff);
  free_zmatrix(U);
  free_zmatrix(V);
  free_zmatrix(a2);
}
