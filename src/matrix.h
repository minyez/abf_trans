#pragma once
#include <cassert>
#include <cstring>
#include <iostream>
#include "linalg.h"

template <typename T>
class matrix
{
public:
    int nr;
    int nc;
    T *c;
    matrix() : nr(0), nc(0), c(nullptr) {}
    matrix(const int &nrows, const int &ncols): nr(nrows), nc(ncols), c(nullptr)
    {
        if (nr&&nc)
            c = new T [nr*nc];
        zero_out();
    }
    matrix(const matrix &m): nr(m.nr), nc(m.nc), c(nullptr)
    {
        if( nr && nc )
        {
            c = new T[nr*nc];
            memcpy(c, m.c, nr*nc*sizeof(T));
        }
    }

    matrix(matrix &&m) : nr(m.nr), nc(m.nc)
    {
        c = m.c;
        m.nr = m.nc = 0;
        m.c = nullptr;
    }

    int size() const { return nc*nr; }
    void zero_out() { for (int i = 0; i < size(); i++) c[i] = 0; }

    T &operator()(const int ir, const int ic) { return c[ir*nc+ic]; }
    const T &operator()(const int ir, const int ic) const { return c[ir*nc+ic]; }

    void operator+=(const matrix<T> &m)
    {
        assert(size() == m.size());
        for (int i = 0; i < size(); i++)
            c[i] += m.c[i];
    }
    void operator+=(const T &cnum)
    {
        for (int i = 0; i < size(); i++)
            c[i] += cnum;
    }

    void operator-=(const matrix<T> &m)
    {
        assert(size() == m.size());
        for (int i = 0; i < size(); i++)
            c[i] -= m.c[i];
    }
    void operator-=(const T &cnum)
    {
        for (int i = 0; i < size(); i++)
            c[i] -= cnum;
    }
    void operator*=(const T &cnum)
    {
        for (int i = 0; i < size(); i++)
            c[i] *= cnum;
    }
    void operator/=(const T &cnum)
    {
        assert(fabs(cnum) > 0);
        for (int i = 0; i < size(); i++)
            c[i] /= cnum;
    }

    void reshape(const int &nrows_new, const int &ncols_new)
    {
        assert ( size() == nrows_new * ncols_new);
        nr = nrows_new;
        nc = ncols_new;
    }

    void resize(const int &nrows_new, const int &ncols_new)
    {
        const int size_new = nrows_new * ncols_new;
        if (size_new)
        {
            if (c)
            {
                if ( size_new != size() )
                {
                    delete [] c;
                    c = new T[size_new];
                }
            }
            else
                c = new T[size_new];
        }
        else
        {
            if(c) delete [] c;
            c = nullptr;
        }
        nr = nrows_new;
        nc = ncols_new;
        zero_out();
    }
    void conj() {};

    void transpose(bool conjugate = false)
    {
        for (int i = 0; i < nr; i++)
            for (int j = i + 1; j < nc; j++)
                c[i*nr + j] = c[j*nc+i];
        int temp = nc;
        nr = temp;
        nc = nr;
        if (conjugate) conj();
    }
};

template <typename T>
void copy(const matrix<T> &src, matrix<T> &dest)
{
    assert(src.size() == dest.size());
    for (int i = 0; i < src.size(); i++)
        dest.c[i] = src.c[i];
}

template <typename T>
matrix<T> operator+(const matrix<T> &m1, const matrix<T> &m2)
{
    assert(m1.nc == m2.nc);
    assert(m1.nr == m2.nr);
    matrix<T> sum = m1;
    sum += m2;
}

template <typename T>
matrix<T> operator+(const matrix<T> &m, const T &cnum)
{
    matrix<T> sum = m;
    sum += cnum;
    return sum;
}

template <typename T>
matrix<T> operator+(const T &cnum, const matrix<T> &m)
{
    return m + cnum;
}

template <typename T>
matrix<T> operator-(const matrix<T> &m1, const matrix<T> &m2)
{
    assert(m1.nc == m2.nc);
    assert(m1.nr == m2.nr);
    matrix<T> sum = m1;
    sum -= m2;
}

template <typename T>
matrix<T> operator-(const matrix<T> &m, const T &cnum)
{
    matrix<T> sum = m;
    sum -= cnum;
    return sum;
}


template <typename T>
matrix<T> operator-(const T &cnum, const matrix<T> &m)
{
    return - m + cnum;
}

template <typename T>
matrix<T> operator*(const matrix<T> &m1, const matrix<T> &m2)
{
    assert(m1.nc == m2.nr);

    matrix<T> prod(m1.nr, m2.nc);
    linalg::gemm('N', 'N', m1.nr, m2.nc, m1.nc,
                 1.0, m1.c, m1.nc, m2.c, m2.nc, 0.0, prod.c, prod.nc);

    return prod;
}

template <typename T>
matrix<T> operator*(const matrix<T> &m, const T &cnum)
{
    matrix<T> sum = m;
    sum *= cnum;
    return sum;
}

template <typename T>
matrix<T> operator*(const T &cnum, const matrix<T> &m)
{
    return m + cnum;
}

template <typename T>
matrix<T> inverse(const matrix<T> &m)
{
    if (m.size() == 0) throw std::invalid_argument("zero size matrix");
    matrix<T> inv;
    inverse(m, inv);
    return inv;
}

template <typename T>
void inverse(const matrix<T> &m, matrix<T> &m_inv)
{
    int lwork = m.nr;
    int info = 0;
    T work[lwork];
    int ipiv[std::min(m.nr, m.nc)];
    m_inv.resize(m.nr, m.nc);
    copy(m, m_inv);
    // debug
    // std::cout << m_inv.size() << " " << m_inv(0, 0) << " " << m_inv(m.nr-1, m.nc-1) << std::endl;
    linalg::getrf(m_inv.nr, m_inv.nc, m_inv.c, m_inv.nc, ipiv, info);
    linalg::getri(m_inv.nr, m_inv.c, m_inv.nc, ipiv, work, lwork, info);
    m_inv.reshape(m_inv.nc, m_inv.nr);
}
