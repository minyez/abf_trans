#pragma once
#include <complex>
#include "base.h"

inline cplxdb get_C_matrix_element(unsigned l, int m, int mp)
{
    if (abs(m) == abs(mp))
    {
        if (m == 0)
            return 1;
        if (m > 0 && mp > 0)
            return std::pow(-1, m)/std::sqrt(2.0);
        if (m > 0 && mp < 0)
            return 1.0/std::sqrt(2.0);
        if (m < 0 && mp > 0)
            return cplxdb(0, -std::pow(-1, mp)/std::sqrt(2.0));
        if (m < 0 && mp < 0)
            return cplxdb(0, 1.0/std::sqrt(2.0));
    }
    return 0;
}

void get_C_matrix(cplxdb *cmat, unsigned l);
