#ifndef GSL_LINALG_H_
#define GSL_LINALG_H_

#include <cmath>

#include "gsl_blas.h"
#include "gsl_matrix.h"
#include "gsl_vector.h"

namespace gsl {

// Non-Block Cholesky.
template <typename T, CBLAS_ORDER O>
void linalg_cholesky_decomp_noblk(matrix<T, O> *A) {
  size_t n = A->size1;
  for (size_t i = 0; i < n; ++i) {
    T l11 = std::sqrt(matrix_get(A, i, i));
    matrix_set(A, i, i, l11);
    if (i + 1 < n) {
      matrix<T, O> l21 = matrix_submatrix(A, i + 1, i, n - i - 1, 1);
      matrix_scale(&l21, 1 / l11);
      matrix<T, O> a22 = matrix_submatrix(A, i + 1, i + 1, n - i - 1,
          n - i - 1);
      blas_syrk(CblasLower, CblasNoTrans, static_cast<T>(-1), &l21,
          static_cast<T>(1), &a22);
    }
  }
}

// Block Cholesky.
//   l11 l11^T = a11
//   l21 = a21 l11^(-T)
//   a22 = a22 - l21 l21^T
//
// Stores result in Lower triangular part.
template <typename T, CBLAS_ORDER O>
void linalg_cholesky_decomp(matrix<T, O> *A) {
  size_t n = A->size1;
  // Block Dimension borrowed from Eigen.
  size_t blk_dim = std::max<size_t>(std::min<size_t>((n / 128) * 16, 8), 128);
  for (size_t i = 0; i < n; i += blk_dim) {
    size_t n11 = std::min<size_t>(blk_dim, n - i);
    matrix<T, O> l11 = matrix_submatrix(A, i, i, n11, n11);
    linalg_cholesky_decomp_noblk(&l11);
    if (i + blk_dim < n) {
      matrix<T, O> l21 = matrix_submatrix(A, i + n11, i, n - i - n11, n11);
      blas_trsm(CblasRight, CblasLower, CblasTrans, CblasNonUnit,
          static_cast<T>(1), &l11, &l21);
      matrix<T, O> a22 = matrix_submatrix(A, i + blk_dim, i + blk_dim, 
          n - i - blk_dim, n - i - blk_dim);
      blas_syrk(CblasLower, CblasNoTrans, static_cast<T>(-1), &l21,
          static_cast<T>(1), &a22);
    }
  }
}

template <typename T, CBLAS_ORDER O>
void linalg_cholesky_svx(const matrix<T, O> *LLT, vector<T> *x) {
  blas_trsv(CblasLower, CblasNoTrans, CblasNonUnit, LLT, x);
  blas_trsv(CblasLower, CblasTrans, CblasNonUnit, LLT, x);
}

}  // namespace gsl

#endif  // GSL_LINALG_H_

