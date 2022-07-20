#pragma once
#include <complex>

using std::complex;

extern "C" {
// BLAS
    float sdot_(const int *N, const float *X, const int *incX, const float *Y,
                const int *incY);
    double ddot_(const int *N, const double *X, const int *incX, const double *Y,
                 const int *incY);
    void cdotc_(complex<float> *result, const int *n, const complex<float> *zx,
                const int *incx, const complex<float> *zy, const int *incy);
    void cdotu_(complex<float> *result, const int *n, const complex<float> *zx,
                const int *incx, const complex<float> *zy, const int *incy);
    void zdotc_(complex<double> *result, const int *n, const complex<double> *zx,
                const int *incx, const complex<double> *zy, const int *incy);
    void zdotu_(complex<double> *result, const int *n, const complex<double> *zx,
                const int *incx, const complex<double> *zy, const int *incy);
    
    void sgemv_(const char *transa, const int *m, const int *n, const float *alpha,
                const float *a, const int *lda, const float *x, const int *incx,
                const float *beta, float *y, const int *incy);
    void dgemv_(const char *transa, const int *m, const int *n, const double *alpha,
                const double *a, const int *lda, const double *x, const int *incx,
                const double *beta, double *y, const int *incy);
    void cgemv_(const char *transa, const int *m, const int *n,
                const complex<float> *alpha, const complex<float> *a,
                const int *lda, const complex<float> *x, const int *incx,
                const complex<float> *beta, complex<float> *y, const int *incy);
    void zgemv_(const char *transa, const int *m, const int *n,
                const complex<double> *alpha, const complex<double> *a,
                const int *lda, const complex<double> *x, const int *incx,
                const complex<double> *beta, complex<double> *y, const int *incy);
    
    void sgemm_(const char *transa, const char *transb, const int *m, const int *n,
                const int *k, const float *alpha, const float *a, const int *lda,
                const float *b, const int *ldb, const float *beta, float *c,
                const int *ldc);
    void dgemm_(const char *transa, const char *transb, const int *m, const int *n,
                const int *k, const double *alpha, const double *a, const int *lda,
                const double *b, const int *ldb, const double *beta, double *c,
                const int *ldc);
    void cgemm_(const char *transa, const char *transb, const int *m, const int *n,
                const int *k, const complex<float> *alpha, const complex<float> *a,
                const int *lda, const complex<float> *b, const int *ldb,
                const complex<float> *beta, complex<float> *c, const int *ldc);
    void zgemm_(const char *transa, const char *transb, const int *m, const int *n,
                const int *k, const complex<double> *alpha,
                const complex<double> *a, const int *lda, const complex<double> *b,
                const int *ldb, const complex<double> *beta, complex<double> *c,
                const int *ldc);
    void dzgemm_(const char *transa, const char *transb, const int *m, const int *n,
                 const int *k, const complex<double> *alpha, const double *a,
                 const int *lda, const complex<double> *b, const int *ldb,
                 const complex<double> *beta, complex<double> *c, const int *ldc);
// LAPACK
    void sgeev_(const char *jobvl, const char *jobvr, const int *n, double *a,
                const int *lda, float *wr, float *wi, float *vl, const int *ldvl,
                float *vr, const int *ldvr, float *work, const int *lwork, int *info);
    void dgeev_(const char *jobvl, const char *jobvr, const int *n, double *a,
                const int *lda, double *wr, double *wi, double *vl, const int *ldvl,
                double *vr, const int *ldvr, double *work, const int *lwork, int *info);
    void cgeev_(const char *jobvl, const char *jobvr, const int *n, float *a,
                const int *lda, complex<float> *wr, complex<float> *wi, complex<float> *vl, const int *ldvl,
                complex<float> *vr, const int *ldvr, complex<float> *work, const int *lwork, int *info);
    void zgeev_(const char *jobvl, const char *jobvr, const int *n, double *a,
                const int *lda, complex<double> *wr, complex<double> *wi, complex<double> *vl, const int *ldvl,
                complex<double> *vr, const int *ldvr, complex<double> *work, const int *lwork, int *info);

    void ssyev_(const char *jobz, const char *uplo, const int *n, float *a,
                const int *lda, float *w, float *work, const int *lwork,
                int *info);
    void dsyev_(const char *jobz, const char *uplo, const int *n, double *a,
                const int *lda, double *w, double *work, const int *lwork,
                int *info);
    void cheev_(const char *jobz, const char *uplo, const int *n,
                complex<float> *a, const int *lda, float *w,
                complex<float> *work, const int *lwork, float *rwork, int *info);
    void zheev_(const char *jobz, const char *uplo, const int *n,
                complex<double> *a, const int *lda, double *w,
                complex<double> *work, const int *lwork, double *rwork, int *info);

    void sgetrf_(const int* m, const int *n, float *A, const int *lda, int *ipiv, int *info);
    void dgetrf_(const int* m, const int *n, double *A, const int *lda, int *ipiv, int *info);
    void cgetrf_(const int* m, const int *n, complex<float> *A, const int *lda, int *ipiv, int *info);
    void zgetrf_(const int* m, const int *n, complex<double> *A, const int *lda, int *ipiv, int *info);

    void sgetri_(const int* n, float *A, const int *lda, int *ipiv, float *work, const int *lwork, int *info);
    void dgetri_(const int* n, double *A, const int *lda, int *ipiv, double *work, const int *lwork, int *info);
    void cgetri_(const int* n, complex<float> *A, const int *lda, int *ipiv, complex<float> *work, const int *lwork, int *info);
    void zgetri_(const int* n, complex<double> *A, const int *lda, int *ipiv, complex<double> *work, const int *lwork, int *info);
}
