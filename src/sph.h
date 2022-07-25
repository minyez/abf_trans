#pragma once
#include "base.h"
#include "matrix.h"

inline cplxdb get_C_matrix_element(int m, int mp)
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

inline cplxdb get_C_matrix_element_aims(int m, int mp)
{
    // TODO: not implemented
    return 0;
}

matrix<cplxdb> get_C_matrix(unsigned l);

matrix<cplxdb> get_C_matrix_aims(unsigned l);
