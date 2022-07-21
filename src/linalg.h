#pragma once
#include "interfaces/blaslapack.h"
#include <stdexcept>

namespace linalg
{
    static inline char revert_trans(const char &trans)
    {
        switch (trans) {
            case 'T':
                return 'N';
            case 't':
                return 'n';
            case 'N':
                return 'T';
            case 'n':
                return 't';
            case 'C':
                throw std::invalid_argument("does not support C, require manually handling");
            case 'c':
                throw std::invalid_argument("does not support C, require manually handling");
            default:
                throw std::invalid_argument("invalid trans character");
        }
    }
    
    static inline char revert_uplo(const char &uplo)
    {
        switch (uplo) {
            case 'U':
                return 'L';
            case 'u':
                return 'l';
            case 'L':
                return 'U';
            case 'l':
                return 'u';
            default:
                throw std::invalid_argument("invalid uplo character");
        }
    }

    template <typename T>
    inline T* transpose(const T* a, const int &n, const int &lda, bool conjugate = false)
    {
        T* a_fort = new T[lda*n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < lda; j++)
                a_fort[i*lda+j] = a[j*n+i];
        return a_fort;
    }

    template <>
    inline complex<float>* transpose(const complex<float>* a, const int &n, const int &lda, bool conjugate)
    {
        complex<float>* a_fort = new complex<float>[lda*n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < lda; j++)
            {
                if (conjugate)
                    a_fort[i*lda+j] = conj(a[j*n+i]);
                else
                    a_fort[i*lda+j] = a[j*n+i];
            }
        return a_fort;
    }

    template <>
    inline complex<double>* transpose(const complex<double>* a, const int &n, const int &lda, bool conjugate)
    {
        complex<double>* a_fort = new complex<double>[lda*n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < lda; j++)
            {
                if (conjugate)
                    a_fort[i*lda+j] = conj(a[j*n+i]);
                else
                    a_fort[i*lda+j] = a[j*n+i];
            }
        return a_fort;
    }

    template <typename T>
    inline void transpose(const T* a_fort, T* a, const int &n, const int &lda, bool conjugate = false)
    {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < lda; j++)
                a[j*n+i] = a_fort[i*lda+j];
    }

    template <>
    inline void transpose(const complex<float>* a_fort, complex<float>* a, const int &n, const int &lda, bool conjugate)
    {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < lda; j++)
                a[j*n+i] = a_fort[i*lda+j];
    }
    
    template <>
    inline void transpose(const complex<double>* a_fort, complex<double>* a, const int &n, const int &lda, bool conjugate)
    {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < lda; j++)
                a[j*n+i] = a_fort[i*lda+j];
    }

    inline float dot(const int &N, const float *X, const int &incX, const float *Y, const int &incY)
    {
        return sdot_(&N, X, &incX, Y, &incY);
    }
    inline double dot(const int &N, const double *X, const int &incX, const double *Y, const int &incY)
    {
        return ddot_(&N, X, &incX, Y, &incY);
    }

    inline complex<float> dot(const int &N, const complex<float> *X, const int &incX,
                              const complex<float> *Y, const int &incY)
    {
        complex<float> res;
        cdotu_(&res, &N, X, &incX, Y, &incY);
        return res;
    }

    inline complex<double> dot(const int &N, const complex<double> *X, const int &incX,
                               const complex<double> *Y, const int &incY)
    {
        complex<double> res;
        zdotu_(&res, &N, X, &incX, Y, &incY);
        return res;
    }
    
    inline complex<float> dotc(const int &N, const complex<float> *X, const int &incX,
                               const complex<float> *Y, const int &incY)
    {
        complex<float> res;
        cdotc_(&res, &N, X, &incX, Y, &incY);
        return res;
    }

    inline complex<double> dotc(const int &N, const complex<double> *X, const int &incX,
                                const complex<double> *Y, const int &incY)
    {
        complex<double> res;
        zdotc_(&res, &N, X, &incX, Y, &incY);
        return res;
    }

    inline void gemm(const char &transa, const char &transb, const int &m, const int &n,
                     const int &k, const float &alpha, const float *a, const int &lda,
                     const float *b, const int &ldb, const float &beta, float *c, const int &ldc)
    {
        sgemm_(&transb, &transa, &n, &m, &k, &alpha, b, &ldb, a, &lda, &beta, c, &ldc);
    }
    
    inline void gemm(const char &transa, const char &transb, const int &m, const int &n,
                     const int &k, const double &alpha, const double *a, const int &lda,
                     const double *b, const int &ldb, const double &beta, double *c, const int &ldc)
    {
        dgemm_(&transb, &transa, &n, &m, &k, &alpha, b, &ldb, a, &lda, &beta, c, &ldc);
    }
    
    inline void gemm(const char &transa, const char &transb, const int &m, const int &n,
                     const int &k, const complex<float> &alpha, const complex<float> *a, const int &lda,
                     const complex<float> *b, const int &ldb, const complex<float> &beta, complex<float> *c, const int &ldc)
    {
        cgemm_(&transb, &transa, &n, &m, &k, &alpha, b, &ldb, a, &lda, &beta, c, &ldc);
    }
    
    inline void gemm(const char &transa, const char &transb, const int &m, const int &n,
                     const int &k, const complex<double> &alpha, const complex<double> *a, const int &lda,
                     const complex<double> *b, const int &ldb, const complex<double> &beta, complex<double> *c, const int &ldc)
    {
        zgemm_(&transb, &transa, &n, &m, &k, &alpha, b, &ldb, a, &lda, &beta, c, &ldc);
    }

    inline void getrf(const int &m, const int &n, float *A, const int &lda, int *ipiv, int &info)
    {
        float *a_fort = transpose(A, n, lda);
        sgetrf_(&m, &n, a_fort, &lda, ipiv, &info);
        transpose(a_fort, A, n, lda);
        delete [] a_fort;
    }
    
    inline void getrf(const int &m, const int &n, double *A, const int &lda, int *ipiv, int &info)
    {
        double *a_fort = transpose(A, n, lda);
        dgetrf_(&m, &n, a_fort, &lda, ipiv, &info);
        transpose(a_fort, A, n, lda);
        delete [] a_fort;
    }
    
    inline void getrf(const int &m, const int &n, complex<float> *A, const int &lda, int *ipiv, int &info)
    {
        complex<float> *a_fort = transpose(A, n, lda);
        cgetrf_(&m, &n, a_fort, &lda, ipiv, &info);
        transpose(a_fort, A, n, lda);
        delete [] a_fort;
    }
    
    inline void getrf(const int &m, const int &n, complex<double> *A, const int &lda, int *ipiv, int &info)
    {
        complex<double> *a_fort = transpose(A, n, lda);
        zgetrf_(&m, &n, a_fort, &lda, ipiv, &info);
        transpose(a_fort, A, n, lda);
        delete [] a_fort;
    }

    inline void getri(const int &n, float *A, const int &lda, int *ipiv, float *work, const int &lwork, int &info)
    {
        float *a_fort = transpose(A, n, lda);
        sgetri_(&n, a_fort, &lda, ipiv, work, &lwork, &info);
        transpose(a_fort, A, n, lda);
        delete [] a_fort;
    }
    
    inline void getri(const int &n, double *A, const int &lda, int *ipiv, double *work, const int &lwork, int &info)
    {
        double *a_fort = transpose(A, n, lda);
        dgetri_(&n, a_fort, &lda, ipiv, work, &lwork, &info);
        transpose(a_fort, A, n, lda);
        delete [] a_fort;
    }
    
    inline void getri(const int &n, complex<float> *A, const int &lda, int *ipiv, complex<float> *work, const int &lwork, int &info)
    {
        complex<float> *a_fort = transpose(A, n, lda);
        cgetri_(&n, a_fort, &lda, ipiv, work, &lwork, &info);
        transpose(a_fort, A, n, lda);
        delete [] a_fort;
    }
    
    inline void getri(const int &n, complex<double> *A, const int &lda, int *ipiv, complex<double> *work, const int &lwork, int &info)
    {
        complex<double> *a_fort = transpose(A, n, lda);
        zgetri_(&n, a_fort, &lda, ipiv, work, &lwork, &info);
        transpose(a_fort, A, n, lda);
        delete [] a_fort;
    }

}
