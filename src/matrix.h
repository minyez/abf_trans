#pragma once
#include <cassert>
#include <vector>
// #include <iostream>
#include <cmath>
#include "base.h"
#include "linalg.h"
#include "vec.h"

template <typename T>
class matrix
{
private:
    int mrank;
public:
    constexpr static const double EQUAL_THRES = DOUBLE_EQUAL_THRES;
    int nr;
    int nc;
    T *c;
    matrix() : nr(0), nc(0), c(nullptr) { mrank = 0; }
    matrix(const int &nrows, const int &ncols): nr(nrows), nc(ncols), c(nullptr)
    {
        if (nr&&nc)
        {
            c = new T [nr*nc];
            mrank = std::min(nr, nc);
        }
        zero_out();
    }
    matrix(const std::vector<vec<T>> &nested_vec): nr(nested_vec.size()), c(nullptr)
    {
        if (nested_vec.size() != 0)
            nc = nested_vec[0].size();
        if (nr&&nc)
        {
            c = new T [nr*nc];
            mrank = std::min(nr, nc);
            for (int ir = 0; ir < nr; ir++)
                for (int ic = 0; ic < std::min(nc, nested_vec[ir].size()); ic++)
                    c[ir*nc+ic] = nested_vec[ir][ic];
        }
    }
    matrix(const int &nrows, const int &ncols, const T * const valarr): nr(nrows), nc(ncols), c(nullptr)
    {
        if (nr&&nc)
        {
            mrank = std::min(nr, nc);
            c = new T [nr*nc];
        }
        zero_out();
        // do not manually check out of bound
        for (int i = 0; i < size(); i++)
            c[i] = valarr[i];
    }
    matrix(const matrix &m): nr(m.nr), nc(m.nc), c(nullptr)
    {
        if( nr && nc )
        {
            c = new T[nr*nc];
            mrank = std::min(nr, nc);
            memcpy(c, m.c, nr*nc*sizeof(T));
        }
    }
    matrix(matrix &&m) : nr(m.nr), nc(m.nc)
    {
        c = m.c;
        mrank = m.mrank;
        m.nr = m.nc = 0;
        m.c = nullptr;
    }

    ~matrix()
    {
        nc = nr = mrank = 0;
        if (c)
        {
            delete [] c;
            c = nullptr;
        }
    }

    int size() const { return nc*nr; }
    void zero_out() { for (int i = 0; i < size(); i++) c[i] = 0; }

    vec<T> get_row(int ir) const
    {
        if (ir < 0 || ir >= nr) throw std::invalid_argument("out-of-bound row index");
        return vec<T>(nc, c+ir*nc);
    }

    void set_diag(const T &v)
    {
        for (int i = 0; i < mrank; i++)
            c[i*nc+i] = v;
    }

    T &operator()(const int ir, const int ic) { return c[ir*nc+ic]; }
    const T &operator()(const int ir, const int ic) const { return c[ir*nc+ic]; }

    matrix<T> & operator=(const matrix<T> &m)
    {
        if (this == &m) return *this;
        resize(m.nr, m.nc);
        memcpy(c, m.c, nr*nc*sizeof(T));
        return *this;
    }

    matrix<T> & operator=(matrix<T> &&m)
    {
        if (this == &m) return *this;
        nr = m.nr;
        nc = m.nc;
        if(c) delete [] c;
        c = m.c;
        m.nr = m.nc = 0;
        m.c = nullptr;
        return *this;
    }

    matrix<T> & operator=(const std::vector<vec<T>> &nested_vec)
    {
        int nr_new = nested_vec.size();
        int nc_new = 0;
        if (nr_new != 0)
            nc_new = nested_vec[0].size();
        resize(nr_new, nc_new);
        if (nr&&nc)
        {
            for (int ir = 0; ir < nr; ir++)
                for (int ic = 0; ic < std::min(nc, nested_vec[ir].size()); ic++)
                    c[ir*nc+ic] = nested_vec[ir][ic];
        }
        return *this;
    }

    matrix<T> & operator=(const T &cnum)
    {
        for (int i = 0; i < size(); i++)
            c[i] = cnum;
        return *this;
    }

    bool operator==(const matrix<T> &m) const
    {
        if (size() == 0 || m.size() == 0) return false;
        if (nc != m.nc && nr != m.nr) return false;
        for (int i = 0; i < size(); i++)
            if (fabs(c[i] - m.c[i]) > matrix<T>::EQUAL_THRES) return false;
        return true;
    }

    void operator+=(const matrix<T> &m)
    {
        assert(size() == m.size());
        for (int i = 0; i < size(); i++)
            c[i] += m.c[i];
    }

    void add_col(const std::vector<T> &v)
    {
        assert(nr == v.size());
        for (int ir = 0; ir < nr; ir++)
            for (int ic = 0; ic < nc; ic++)
                c[ir*nc+ic] += v[ir];
    }

    void operator+=(const std::vector<T> &v)
    {
        assert(nc == v.size());
        for (int i = 0; i < nr; i++)
            for (int ic = 0; ic < nc; ic++)
                c[i*nc+ic] += v[ic];
    }

    void operator+=(const vec<T> &v)
    {
        assert(nc == v.size());
        for (int i = 0; i < nr; i++)
            for (int ic = 1; ic < nc; ic++)
                c[i*nc+ic] += v.c[ic];
    }

    void operator-=(const std::vector<T> &v)
    {
        assert(nc == v.size());
        for (int i = 0; i < nr; i++)
            for (int ic = 1; ic < nc; ic++)
                c[i*nc+ic] -= v[ic];
    }

    void operator-=(const vec<T> &v)
    {
        assert(nc == v.size());
        for (int i = 0; i < nr; i++)
            for (int ic = 1; ic < nc; ic++)
                c[i*nc+ic] -= v.c[ic];
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
        mrank = std::min(nr, nc);
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
        mrank = std::min(nr, nc);
        zero_out();
    }
    void conj() {};

    void transpose(bool conjugate = false)
    {
        for (int i = 0; i < nr; i++)
            for (int j = i + 1; j < nc; j++)
            {
                T temp = c[i*nr + j];
                c[i*nr + j] = c[j*nc+i];
                c[j*nc+i] = temp;
            }
        int temp = nc;
        nr = temp;
        nc = nr;
        if (conjugate) conj();
    }

    T det() const { return get_determinant(*this); }
};

template <> inline void matrix<complex<float>>::conj()
{
    for (int i = 0; i < size(); i++)
        c[i] = std::conj(c[i]);
}

template <> inline void matrix<complex<double>>::conj()
{
    for (int i = 0; i < size(); i++)
        c[i] = std::conj(c[i]);
}

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
    return sum;
}

template <typename T>
matrix<T> operator+(const matrix<T> &m, const std::vector<T> &v)
{
    assert(m.nc == v.size());
    matrix<T> mnew(m);
    mnew += v;
    return mnew;
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
    assert(m1.nc == m2.nc && m1.nr == m2.nr);
    matrix<T> mnew = m1;
    mnew -= m2;
    return mnew;
}

template <typename T>
matrix<T> operator-(const matrix<T> &m, const std::vector<T> &v)
{
    assert(m.nc == v.size());
    matrix<T> mnew = m;
    mnew -= v;
    return mnew;
}

template <typename T>
matrix<T> operator-(const matrix<T> &m, const T &cnum)
{
    matrix<T> mnew = m;
    mnew -= cnum;
    return mnew;
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
vec<T> operator*(const matrix<T> &m, const vec<T> &v)
{
    assert(m.nc == v.n);
    vec<T> mv(m.nr);
    linalg::gemv('N', m.nr, m.nc, 1.0, m.c, m.nc, v.c, 1, 0.0, mv.c, 1);
    return mv;
}

template <typename T>
vec<T> operator*(const vec<T> &v, const matrix<T> &m)
{
    assert(m.nr == v.n);
    vec<T> mv(m.nc);
    /* linalg::gemv('N', ); */
    linalg::gemv('T', m.nr, m.nc, 1.0, m.c, m.nc, v.c, 1, 0.0, mv.c, 1);
    return mv;
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

template <typename T>
matrix<T> transpose(const matrix<T> &m, bool conjugate = false)
{
    matrix<T> mnew(m);
    mnew.transpose(conjugate);
    return mnew;
}

template <typename T>
std::ostream & operator<<(std::ostream & os, const matrix<T> &m)
{
    for (int i = 0; i < m.nr; i++)
    {
        for (int j = 0; j < m.nc; j++)
            os << m(i, j) << " ";
        os << std::endl;
    }
    return os;
}

template <typename T>
matrix<complex<T>> to_complex(const matrix<T> &m)
{
    matrix<complex<T>> cm(m.nr, m.nc);
    for (int i = 0; i < m.size(); i++)
        cm.c[i] = m.c[i];
    return cm;
}

template <typename T>
matrix<T> get_real(const matrix<complex<T>> &cm)
{
    matrix<T> m(cm.nr, cm.nc);
    for (int i = 0; i < cm.size(); i++)
        m.c[i] = cm.c[i].real();
    return m;
}

template <typename T>
matrix<T> get_imag(const matrix<complex<T>> &cm)
{
    matrix<T> m(cm.nr, cm.nc);
    for (int i = 0; i < cm.size(); i++)
        m.c[i] = cm.c[i].imag();
    return m;
}

template <typename T>
T get_determinant(matrix<T> m)
{
    int lwork = m.nr;
    int info = 0;
    int mrank = std::min(m.nr, m.nc);
    T work[lwork];
    int ipiv[mrank];
    // debug
    linalg::getrf(m.nr, m.nc, m.c, m.nc, ipiv, info);
    T det = 1;
    for (int i = 0; i < mrank; i++)
    {
        /* std::cout << i << " " << ipiv[i] << " " << m.c[i*m.nc+i] << " "; */
        det *= (2*int(ipiv[i] == (i+1))-1) * m.c[i*m.nc+i];
    }
    return det;
}

template <typename T>
matrix<double> to_double(const matrix<T> &mat)
{
    matrix<double> dmat(mat.nr, mat.nc);
    if (dmat.size())
    {
        for (int ir = 0; ir < mat.nr; ir++)
            for (int ic = 0; ic < mat.nc; ic++)
                dmat(ir, ic) = mat(ir, ic);
    }
    return dmat;
}
