#pragma once
#include "base.h"
#include "matrix.h"
#include <string>
#include <cassert>

// real spherical harmonics used in the original Blanco's paper
inline cplxdb get_C_matrix_element_orig(int m, int mp)
{
    if (abs(m) == abs(mp))
    {
        if (m == 0)
            return cplxdb(1.0, 0.0);
        if (m > 0 && mp > 0)
            return cplxdb(std::pow(-1.0, m)/std::sqrt(2.0), 0.0);
        if (m > 0 && mp < 0)
            return cplxdb(1.0/std::sqrt(2.0), 0.0);
        if (m < 0 && mp > 0)
            return cplxdb(0, -std::pow(-1, mp)/std::sqrt(2.0));
        if (m < 0 && mp < 0)
            return cplxdb(0, 1.0/std::sqrt(2.0));
    }
    return cplxdb(0., 0.);
}

// real spherical harmonics used in the aims code
inline cplxdb get_C_matrix_element_aims(int m, int mp)
{
    if (abs(m) == abs(mp))
    {
        if (m == 0)
            return cplxdb(1.0, 0.0);
        if (m > 0 && mp > 0)
            return cplxdb(1.0/std::sqrt(2.0), 0.0);
        if (m > 0 && mp < 0)
            return cplxdb(std::pow(-1.0, m)/std::sqrt(2.0), 0.0);
        if (m < 0 && mp > 0)
            return cplxdb(0, std::pow(-1, mp)/std::sqrt(2.0));
        if (m < 0 && mp < 0)
            return cplxdb(0, -1.0/std::sqrt(2.0));
    }
    return cplxdb(0., 0.);
}

inline cplxdb get_C_matrix_element(int m, int mp, const std::string &choice = "orig")
{
    if (choice == "orig")
        return get_C_matrix_element_orig(m, mp);
    if (choice == "aims")
        return get_C_matrix_element_aims(m, mp);
    else
        throw std::invalid_argument("Unknown choice of real spherical harmonics");
}

matrix<cplxdb> get_C_matrix(unsigned l, const std::string &choice = "orig");
