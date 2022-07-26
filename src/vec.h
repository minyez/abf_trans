#include "base.h"
#include <cstring>
#include <complex>
#include <cassert>
#include <cmath>

template <typename T>
class vec
{
public:
    int n;
    T* c;
    constexpr static const double EQUAL_THRES = DOUBLE_EQUAL_THRES;

    vec(): n(0), c(nullptr) {};
    vec(const int &n_in): n(n_in) { if (n_in > 0) c = new T [n]; }
    vec(const int &n_in, const T * const valarr): n(n_in), c(nullptr)
    {
        if (n>0)
        {
            c = new T [n];
        }
        zero_out();
        // do not manually check out of bound
        for (int i = 0; i < size(); i++)
            c[i] = valarr[i];
    }
    vec(const vec &v): n(v.n), c(nullptr)
    {
        if(n > 0)
        {
            c = new T[n];
            memcpy(c, v.c, n*sizeof(T));
        }
    }
    vec(vec &&v) : n(v.n)
    {
        c = v.c;
        v.n = 0;
        v.c = nullptr;
    }
    ~vec() { n = 0; if (c) { delete [] c; c = nullptr;}}

    int size() const { return n; }

    void zero_out() { for (int i = 0; i < size(); i++) c[i] = 0; }

    T &operator[](const int i) { return c[i]; }
    const T &operator[](const int i) const { return c[i]; }

    vec<T> & operator=(const vec<T> &v)
    {
        if (this == &v) return *this;
        resize(v.n);
        memcpy(c, v.c, n*sizeof(T));
        return *this;
    }

    vec<T> & operator=(vec<T> &&v)
    {
        if (this == &v) return *this;
        n = v.n;
        if(c) delete [] c;
        c = v.c;
        v.c = nullptr;
        return *this;
    }

    bool operator<(const vec<T> &v) const
    {
        return norm2(*this) < norm2(v);
    }

    bool operator>(const vec<T> &v) const
    {
        return norm2(*this) > norm2(v);
    }

    void resize(const int &n_new)
    {
        const int size_new = n_new;
        if (size_new > 0)
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
        n = n_new;
        zero_out();
    }

    bool operator==(const vec<T> &v) const
    {
        if (size() == 0 || v.size() == 0) return false;
        for (int i = 0; i < size(); i++)
            if (fabs(c[i] - v.c[i]) > vec<T>::EQUAL_THRES) return false;
        return true;
    }

    void operator+=(const T &cnum)
    {
        for (int i = 0; i < size(); i++)
            c[i] += cnum;
    }

    void operator+=(const vec<T> &v)
    {
        assert(n == v.n);
        for (int i = 0; i < size(); i++)
            c[i] += v.c[i];
    }

    void operator-=(const T &cnum)
    {
        for (int i = 0; i < size(); i++)
            c[i] -= cnum;
    }

    void operator-=(const vec<T> &v)
    {
        assert(n == v.n);
        for (int i = 0; i < size(); i++)
            c[i] -= v.c[i];
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
};

template <typename T>
vec<T> conj(const vec<T> &v)
{
    vec<T> conjv(v);
    return conjv;
}

template <> inline vec<std::complex<float>> conj(const vec<std::complex<float>> &v)
{
    auto conjv(v);
    for (int i = 0; i < conjv.size(); i++)
        conjv.c[i] = std::conj(conjv.c[i]);
    return conjv;
}

template <> inline vec<std::complex<double>> conj(const vec<std::complex<double>> &v)
{
    auto conjv(v);
    for (int i = 0; i < conjv.size(); i++)
        conjv.c[i] = std::conj(conjv.c[i]);
    return conjv;
}

template <typename T>
vec<T> operator+(const vec<T> &v, const T &cnum)
{
    vec<T> sum = v;
    sum += cnum;
    return sum;
}

template <typename T>
vec<T> operator+(const T &cnum, const vec<T> &v)
{
    return v + cnum;
}

template <typename T>
vec<T> operator+(const vec<T> &v1, const vec<T> &v2)
{
    assert(v1.n == v2.n);
    vec<T> sum = v1;
    sum += v2;
    return sum;
}

template <typename T>
vec<T> operator-(const vec<T> &v1, const vec<T> &v2)
{
    assert(v1.n == v2.n);
    vec<T> vnew = v1;
    vnew -= v2;
    return vnew;
}

template <typename T>
vec<T> operator-(const vec<T> &v, const T &cnum)
{
    vec<T> mnew = v;
    mnew -= cnum;
    return mnew;
}


template <typename T>
vec<T> operator-(const T &cnum, const vec<T> &v)
{
    return - v + cnum;
}

template <typename T>
std::ostream & operator<<(std::ostream & os, const vec<T> &v)
{
    for (int i = 0; i < v.n; i++)
        os << v[i] << " ";
    os << std::endl;
    return os;
}

template <typename T>
vec<T> operator*(const vec<T> &v, const T &cnum)
{
    vec<T> mnew = v;
    mnew *= cnum;
    return mnew;
}


template <typename T>
vec<T> operator*(const T &cnum, const vec<T> &v)
{
    return v * cnum;
}

template <typename T>
T norm2(vec<T> v)
{
    return norm(v.c, v.n, 2);
}

template <typename T>
T norm2(vec<std::complex<T>> v)
{
    return norm(v.c, v.n, 2);
}
