#define _USE_MATH_DEFINES
#include "cmatrix.h"
#include "dsyev.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  int n = 5;
  int loop = 4; // a^(2^loop) を計算

  double **a, **a2, **a4; /* matrix */
  double *v, *v2, *tmp;
  a = alloc_dmatrix(n, n);
  a2 = alloc_dmatrix(n, n);
  a4 = alloc_dmatrix(n, n);
  v = alloc_dvector(n);
  v2 = alloc_dvector(n);
  tmp = alloc_dvector(n);

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      a[i][j] = 0.0;
  for (int i = 0; i < n - 1; i++) {
    a[i][i] = 2;
    a[i][i + 1] = a[i + 1][i] = -1;
  }
  a[n - 1][n - 1] = 1;
  printf("Matrix A:\n");
  fprint_dmatrix(stdout, n, n, a);

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      a2[i][j] = a[i][j];
  for (int l = 0; l < loop; l++) {
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) {
        a4[i][j] = 0.0;
        for (int k = 0; k < n; k++) {
          a4[i][j] += a2[i][k] * a2[k][j];
        }
      }
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        a2[i][j] = a4[i][j];
  }
  printf("Matrix A^(2^%d):\n", loop);
  fprint_dmatrix(stdout, n, n, a4);
  for (int i = 0; i < n; i++) {
    v[i] = 1.0;
  }
  for (int i = 0; i < n; i++) {
    tmp[i] = 0.0;
    for (int j = 0; j < n; j++)
      tmp[i] += a4[i][j] * v[j];
  }
  for (int i = 0; i < n; i++)
    v[i] = tmp[i];
  for (int i = 0; i < n; i++) {
    v2[i] = 0.0;
    for (int j = 0; j < n; j++)
      v2[i] += a[i][j] * v[j];
  }
  double d = 0.0, m = 0.0;
  for (int i = 0; i < n; i++) {
    d += v[i] * v2[i];
    m += v2[i] * v2[i];
  }
  double ans1 = m / d;
  double ans2 = 2.0 * (1.0 - cos((double)(2 * n - 1) * M_PI / (2 * n + 1)));
  printf("\nlambda max: %10.10f\n", ans1);
  printf("lambda max(theory): %10.10f\n", ans2);
  printf("diff:%.10e", sqrt((ans1 - ans2) * (ans1 - ans2)));

  free_dmatrix(a);
  free_dmatrix(a2);
  free_dmatrix(a4);
  free_dvector(v);
  free_dvector(v2);
  free_dvector(tmp);
}
