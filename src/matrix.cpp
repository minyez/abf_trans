#include "matrix.h"
#include <complex>
#include <cassert>
#include <cmath>
using std::complex;

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

