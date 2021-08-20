#include "cmatrix.h"
#include "dsyev.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  char *filename;
  FILE *fp;

  int m, n;
  double **a, **a2; /* matrix */
  double *w;        /* eigenvalues */
  double *work;     /* working area */

  int lwork, info;
  char jobz = 'V'; /* Compute eigenvalues and eigenvectors */
  char uplo = 'U'; /* Upper triangle of A is stored */

  if (argc < 2) {
    fprintf(stderr, "Usage: %s inputfile\n", argv[0]);
    exit(1);
  }
  filename = argv[1];

  /* read matrix A */
  fp = fopen(filename, "r");
  if (fp == NULL) {
    fprintf(stderr, "Error: file can not open\n");
    exit(1);
  }
  read_dmatrix(fp, &m, &n, &a);
  if (m != n) {
    fprintf(stderr, "Error: matrix is not square\n");
    exit(1);
  }
  printf("Matrix A:\n");
  fprint_dmatrix(stdout, n, n, a);

  a2 = alloc_dmatrix(n, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      a2[i][j] = a[i][j];

  /* perform eigenvalue decomposition */
  w = alloc_dvector(n);
  lwork = 3 * n - 1;
  work = alloc_dvector(lwork);
  dsyev_(&jobz, &uplo, &n, mat_ptr(a), &n, vec_ptr(w), vec_ptr(work), &lwork, &info);
  if (info != 0) {
    fprintf(stderr, "Error: LAPACK::dsyev failed\n");
    exit(1);
  }
  printf("Eigenvalues:\n");
  fprint_dvector(stdout, n, w);
  printf("Eigenvectors [each column represents each eigenvector]:\n");
  fprint_dmatrix(stdout, n, n, a);

  double **U;
  U = alloc_dmatrix(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      U[i][j] = 0.0;
      for (int k = 0; k < n; k++) {
        U[i][j] += a[i][k] * a[j][k];
      }
      if (j == i) U[i][j] -= 1.0;
    }
  }
  printf("Matrix (tU)U-I:\n");
  fprint_dmatrix_e(stdout, n, n, U);

  double **V;
  V = alloc_dmatrix(n, n);

  // for (int i = 0; i < n; i++) {
  //   double diff = 0.0;
  //   for (int j = 0; j < n; j++) {
  //     V[j][0] = 0.0;
  //     for (int k = 0; k < n; k++) {
  //       V[j][0] += a2[j][k] * a[k][i];
  //     }
  //     V[j][0] -= w[i] * a[j][i];
  //     diff += V[j][0] * V[j][0];
  //   }
  //   printf("%d : diff=%f\n", i, sqrt(diff));
  // }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      U[i][j] = 0.0;
      for (int k = 0; k < n; k++) {
        U[i][j] += a[i][k] * a2[k][j];
      }
    }
  }
  printf("Matrix V:\n");
  fprint_dmatrix(stdout, n, n, U);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      // V[i][j] = 0.0;
      // for (int k = 0; k < n; k++) {
      //   V[i][j] += U[i][k] * U[j][k];
      // }
      V[i][j] = a[i][j] * w[i];
    }
  }
  printf("Matrix V':\n");
  fprint_dmatrix(stdout, n, n, V);
  // for (int i = 0; i < n; i++) {
  //   w[i] *= w[i];
  // }
  // fprint_dvector(stdout, n, w);

  double *diff;
  diff = alloc_dvector(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      diff[i] += (U[i][j] - V[i][j]) * (U[i][j] - V[i][j]);
    }
    diff[i] = sqrt(diff[i]);
  }
  printf("diffs:\n");
  fprint_dvector(stdout, n, diff);

  free_dmatrix(a);
  free_dvector(w);
  free_dvector(diff);
  free_dmatrix(U);
  free_dmatrix(V);
  // free_dmatrix(D);
  free_dmatrix(a2);
}
