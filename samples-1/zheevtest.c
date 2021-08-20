/*
 * SIZE*SIZE型のエルミート行列の固有値と固有ベクトルを計算
 *  (1  -2i)
 *  (+2i 1)
 */

#include <stdio.h>
#include <complex.h>
#define SIZE 2

int main(void) {
// 複素行列(エルミート)の対角化 (入力行列の配列は対角化後にユニタリ行列になる)
  char   jobz = 'V'           ;// 固有ベクトルを計算する
  char   uplo = 'U'           ;// Aに上三角行列を格納
  int    n    = SIZE          ;// 対角化する正方行列のサイズ

  double _Complex A[SIZE*SIZE];// 対角化する行列。対角化後は固有ベクトルが並ぶ
  A[0]=1; A[2]=-2*I;
  A[1]=2*I; A[3]=1;

  double w[SIZE]              ;// 実固有値が小さい順に入る

  int    lda  = SIZE          ;// 対角化する正方行列のサイズ
  double _Complex work[6*SIZE];// 対角化する際に使用するメモリ
  int    lwork = 6*SIZE       ;// workの次元
  double rwork[3*SIZE-2]      ;// 3*SIZE-2で固定
  int    info                 ;// 成功すれば0、失敗すれば0以外を返す

  zheev_( &jobz, &uplo, &n, A, &lda, w, work, &lwork, rwork, &info );
  // 1番目の固有値 : w[0]  1番目の固有ベクトル : (A[0] A[1])
  // 2番目の固有値 : w[1]  2番目の固有ベクトル : (A[2] A[3])

  printf("first eigenvalue:%5.3lf\n",w[0]);
  printf("first eigenvector:(%5.3lf %+5.3lf*I %5.3lf %+5.3lf*I)\n",
      creal(A[0]), cimag(A[0]), creal(A[1]), cimag(A[1]));

  printf("second eigenvalue:%5.3lf\n",w[1]);
  printf("second eigenvector:(%5.3lf %+5.3lf*I %5.3lf %+5.3lf*I)\n",
      creal(A[2]), cimag(A[2]), creal(A[3]), cimag(A[3]));

  return 0;
}