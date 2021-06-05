#include "samples-1/cmatrix.h"
#include "samples-1/mersenne_twister.h"
#include <complex.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define I _Complex_I
#define THRESHOLD 1e-12 /* 直交行列、ユニタリ行列の判定に用いる閾値 */
double range = 100.0;

/* ランダムな n 次実対称行列を作る。各要素 ∈ (-100,100)。*/
void genSymMatRand(int n, double **mat) {
  for (int i = 0; i < n; i++)
    for (int j = 0; j <= i; j++)
      mat[i][j] = mat[j][i] = genrand_real3() * 2.0 * range - range;
}

/* n次正方行列 mat が対称行列なら 1、そうでないなら 0 を返す */
int isSymMat(int n, const double **mat) {
  int flg = 1;
  for (int i = 0; i < n; i++)
    for (int j = 0; j <= i; j++)
      if (mat[i][j] != mat[j][i])
        flg = 0;

  return flg;
}

/* ランダムな n 次エルミート行列を作る。各要素の実部、虚部 ∈ (-100,100)。*/
void genHrmMatRand(int n, double complex **zmat) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < i; j++) {
      zmat[i][j] = (genrand_real3() * 2.0 * range - range) + (genrand_real3() * 2.0 * range - range) * I;
      zmat[j][i] = conj(zmat[i][j]);
    }
    zmat[i][i] = genrand_real3() * 2.0 * range - range;
  }
}

/* n次正方行列 zmat がエルミートなら 1、そうでないなら 0 を返す */
int isHrmMat(int n, const double complex **zmat) {
  int flg = 1;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      if (zmat[i][j] != conj(zmat[j][i]))
        flg = 0;
    }
  }
  return flg;
}

/*
 * n 個の順列をあらわす perm に次の順列を入れる。
 * 渡されたのが最後の順列なら 0 を、次があれば 1 を返す。
 * sgn には置換の符号を入れる。
 */
int nextPerm(int n, int *perm, int *sgn) {
  int a = n - 2;
  while (0 <= a && perm[a + 1] <= perm[a])
    --a;
  if (a < 0)
    return 0;
  int b = n - 1;
  while (perm[b] <= perm[a])
    --b;
  int t = perm[b];
  perm[b] = perm[a], perm[a] = t;
  *sgn *= -1;
  for (int i = a + 1, j = n - 1; i < j; ++i, --j) {
    t = perm[i], perm[i] = perm[j], perm[j] = t;
    *sgn *= -1;
  }
  return 1;
}

/* 実行列行列式 */
double det(int n, const double **mat) {
  double sum = 0;
  int *perm = alloc_ivector(n);
  for (int i = 0; i < n; i++)
    perm[i] = i;
  int sgn = 1;
  do {
    double prod = sgn;
    for (int i = 0; i < n; i++)
      prod *= mat[i][perm[i]];
    sum += prod;
  } while (nextPerm(n, perm, &sgn));
  free_ivector(perm);
  return sum;
}

/* 複素行列行列式 */
double complex zdet(int n, const double complex **zmat) {
  double complex sum = 0;
  int *perm = alloc_ivector(n);
  for (int i = 0; i < n; i++)
    perm[i] = i;
  int sgn = 1;
  do {
    double complex prod = sgn;
    for (int i = 0; i < n; i++)
      prod *= zmat[i][perm[i]];
    sum += prod;
  } while (nextPerm(n, perm, &sgn));
  free_ivector(perm);
  return sum;
}

/* m × n 実行列の 1 ノルムを計算。 */
double norm1(int m, int n, const double **mat) {
  double ret = 0;
  for (int j = 0; j < n; j++) {
    double sum = 0;
    for (int i = 0; i < m; i++)
      sum += fabs(mat[i][j]);
    ret = fmax(ret, sum);
  }
  return ret;
}

/* m × n 実行列の ∞ ノルムを計算。 */
double normInf(int m, int n, const double **mat) {
  double ret = 0;
  for (int i = 0; i < m; i++) {
    double sum = 0;
    for (int j = 0; j < n; j++)
      sum += fabs(mat[i][j]);
    ret = fmax(ret, sum);
  }
  return ret;
}

/* m × n 複素行列の 1 ノルムを計算。 */
double znorm1(int m, int n, const double complex **mat) {
  double ret = 0.0;
  for (int j = 0; j < n; j++) {
    double sum = 0.0;
    for (int i = 0; i < m; i++)
      sum += cabs(mat[i][j]);
    ret = fmax(ret, sum);
  }
  return ret;
}

/* m × n 複素行列の ∞ ノルムを計算。 */
double znormInf(int m, int n, const double complex **mat) {
  double ret = 0.0;
  for (int i = 0; i < m; i++) {
    double sum = 0.0;
    for (int j = 0; j < n; j++)
      sum += cabs(mat[i][j]);
    ret = fmax(ret, sum);
  }
  return ret;
}

/* 実ベクトルの内積 */
double inprod(int n, const double *u, const double *v) {
  double ret = 0.0;
  for (int i = 0; i < n; i++)
    ret += u[i] * v[i];
  return ret;
}

/* 複素ベクトルの内積 */
double complex zinprod(int n, const double complex *u, const double complex *v) {
  double complex ret = 0.0;
  for (int i = 0; i < n; i++)
    ret += conj(u[i]) * v[i];
  return ret;
}

/* ランダムな n 次直交行列を作る */
void genOrthoMatRand(int n, double **mat) {
  /* とりあえずランダムな n 次ベクトルを n 個作る。各要素 ∈ (-100,100)。
     n 個が線形独立でないなら作り直す */
  do {
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        mat[i][j] = genrand_real3() * 2.0 * range - range;
    printf("trial\n");
  } while (fabs(det(n, mat)) < THRESHOLD);

  /* Gram-Schmidt の直交化法 */
  double *u2 = alloc_dvector(n);  /* よく使うので既にできた基底の自身との内積を保持 */
  double *tmp = alloc_dvector(n); /* 引き算する値の保持。まとめて引く。 */
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      tmp[j] = 0;
    for (int j = 0; j < i; j++) {
      double c = inprod(n, mat[i], mat[j]) / u2[j];
      for (int k = 0; k < n; k++)
        tmp[k] += c * mat[j][k];
    }
    for (int j = 0; j < n; j++)
      mat[i][j] -= tmp[j];
    u2[i] = inprod(n, mat[i], mat[i]);
  }
  /* 正規化 */
  for (int i = 0; i < n; i++) {
    u2[i] = sqrt(u2[i]);
    for (int j = 0; j < n; j++)
      mat[i][j] /= u2[i];
  }
  free_dvector(u2);
  free_dvector(tmp);
}

/* n 次正方行列 mat が直交行列なら 1、そうでないなら 0 を返す */
int isOrthoMat(int n, const double **mat) {
  double **prod = alloc_dmatrix(n, n);

  /* mat と mat の転置の行列積を計算 */
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      double sum = 0;
      for (int k = 0; k < n; k++)
        sum += mat[i][k] * mat[j][k];
      prod[i][j] = sum;
    }

  /* n 次単位行列との差の行列のノルムを計算 */
  for (int i = 0; i < n; i++)
    prod[i][i] -= 1.0;
  double nrm = norm1(n, n, prod);

  printf("1norm diff = %.30lf\n", nrm);
  printf("inf norm diff = %.30lf\n", normInf(n, n, prod));
  free_dmatrix(prod);
  /* 単位行列との差のノルムが閾値未満なら直交行列と認める */
  if (nrm < THRESHOLD) return 1;
  return 0;
}

/* ランダムな n 次ユニタリ行列を作る */
void genUnitMatRand(int n, double complex **mat) {
  /* とりあえずランダムな n 次複素ベクトルを n 個作る。各要素 ∈ (-100,100)。
     n 個が線形独立でないなら作り直す */
  do {
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        mat[i][j] = (genrand_real3() * 2.0 * range - range) + (genrand_real3() * 2.0 * range - range) * I;
    printf("trial\n");
  } while (cabs(zdet(n, mat)) < THRESHOLD);

  /* Gram-Schmidt の直交化法 */
  double *u2 = alloc_dvector(n);          /* よく使うので既にできた基底の自身との内積を保持 */
  double complex *tmp = alloc_zvector(n); /* 引き算する値の保持。まとめて引く。 */
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      tmp[j] = 0;
    for (int j = 0; j < i; j++) {
      double complex c = zinprod(n, mat[j], mat[i]) / u2[j];
      for (int k = 0; k < n; k++)
        tmp[k] += c * mat[j][k];
    }
    for (int j = 0; j < n; j++)
      mat[i][j] -= tmp[j];
    /* 複素ベクトルでも自身との内積は実数のはず */
    u2[i] = creal(zinprod(n, mat[i], mat[i]));
  }
  /* 正規化 */
  for (int i = 0; i < n; i++) {
    u2[i] = sqrt(u2[i]);
    for (int j = 0; j < n; j++)
      mat[i][j] /= u2[i];
  }
  free_dvector(u2);
  free_zvector(tmp);
}

/* n 次正方行列 mat がユニタリ行列なら 1、そうでないなら 0 を返す */
int isUnitMat(int n, const double complex **mat) {
  double complex **prod = alloc_zmatrix(n, n);

  /* mat と mat のエルミート共役の行列積を計算 */
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      double complex sum = 0.0;
      for (int k = 0; k < n; k++)
        sum += conj(mat[k][i]) * mat[k][j];
      prod[i][j] = sum;
    }

  /* n 次単位行列との差の行列のノルムを計算 */
  for (int i = 0; i < n; i++)
    prod[i][i] -= 1.0;
  double nrm = znorm1(n, n, prod);

  printf("1norm diff = %.30lf\n", nrm);
  printf("inf norm diff = %.30lf\n", znormInf(n, n, prod));
  free_zmatrix(prod);
  /* 単位行列との差のノルムが閾値未満ならユニタリ行列と認める */
  if (nrm < THRESHOLD) return 1;
  return 0;
}

int main(int argc, char *argv[]) {
  int n = 5;
  /* 1つ目のコマンドライン引数は行列のサイズ n（デフォルトは 5） */
  if (argc >= 2) n = atoi(argv[1]);
  /* 2つ目のコマンドラインはランダムな行列要素の範囲 (-range,range)（デフォルトは (-100,100） */
  if (argc >= 3) range = atof(argv[2]);
  init_genrand((unsigned long)time(NULL));

  double **mat;
  mat = alloc_dmatrix(n, n);
  double complex **zmat;
  zmat = alloc_zmatrix(n, n);

  /* 実対称行列 */
  genSymMatRand(n, mat);
  printf("mat = ");
  fprint_dmatrix(stdout, n, n, mat);
  printf(isSymMat(n, mat) ? "Symmetry matrix\n" : "Non-Symmetry matrix\n");
  printf("\n");

  /* エルミート行列 */
  genHrmMatRand(n, zmat);
  printf("zmat = ");
  fprint_zmatrix(stdout, n, n, zmat);
  printf(isHrmMat(n, zmat) ? "Hermitian matrix\n" : "Non-Hermitian matrix\n");
  printf("\n");

  /* 直交行列 */
  genOrthoMatRand(n, mat);
  printf("mat = ");
  fprint_dmatrix(stdout, n, n, mat);
  printf(isOrthoMat(n, mat) ? "Orthogonal matrix\n" : "Non-Orthogonal matrix\n");
  printf("\n");

  /* ユニタリ行列 */
  genUnitMatRand(n, zmat);
  printf("zmat = ");
  fprint_zmatrix(stdout, n, n, zmat);
  printf(isUnitMat(n, zmat) ? "Unitary matrix\n" : "Non-Unitary matrix\n");

  free_dmatrix(mat);
  free_zmatrix(zmat);
  return 0;
}