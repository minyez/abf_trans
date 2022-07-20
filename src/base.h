#pragma once
#include <complex>

typedef std::complex<double> cplxdb;

inline size_t get_msize(unsigned l)
{
    return 2*l+1;
}

