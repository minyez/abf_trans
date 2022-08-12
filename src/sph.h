#pragma once
#include "base.h"
#include "matrix.h"
#include <string>
#include <cassert>

// real spherical harmonics used in the original Blanco's paper
inline cplxdb get_C_matrix_element_orig(int m, int mp)
{
    const double sqrt2 = M_SQRT2;
    if (abs(m) == abs(mp))
    {
        if (m == 0)
            return cplxdb(1.0, 0.0);
        if (m > 0 && mp > 0)
            return cplxdb(std::pow(-1.0, m)/sqrt2, 0.0);
        if (m > 0 && mp < 0)
            return cplxdb(1.0/sqrt2, 0.0);
        if (m < 0 && mp > 0)
            return cplxdb(0, -std::pow(-1, mp)/sqrt2);
        if (m < 0 && mp < 0)
            return cplxdb(0, 1.0/sqrt2);
    }
    return cplxdb(0., 0.);
}

// real spherical harmonics used in the aims code
inline cplxdb get_C_matrix_element_aims(int m, int mp)
{
    const double sqrt2 = M_SQRT2;
    if (abs(m) == abs(mp))
    {
        if (m == 0)
            return cplxdb(1.0, 0.0);
        if (m > 0 && mp > 0)
            return cplxdb(1.0/sqrt2, 0.0);
        if (m > 0 && mp < 0)
            return cplxdb(std::pow(-1.0, m)/sqrt2, 0.0);
        if (m < 0 && mp > 0)
            return cplxdb(0, -std::pow(-1, mp)/sqrt2);
        if (m < 0 && mp < 0)
            return cplxdb(0, 1.0/sqrt2);

        // original code from src/sym_base.f90
        // if ((abs(m).ne.m).and. (abs(m2) .ne. m2))then
        //     ! upper left, m<0 + m2<0
        //     element = dcmplx(0.d0,1.d0)
        // elseif ((abs(m) == m) .and. (abs(m2) .ne. m2))then
        //     ! lower left, m>0 + m2<0
        //     element = dcmplx(1.d0,0.d0)
        // elseif ((abs(m) == m) .and. (abs(m2) == m2))then
        //     ! lower right, m>0 + m2>0
        //     element = dcmplx((-1)**(l-(l-m)),0.d0)
        // elseif (.not.(abs(m).eq.m) .and. (abs(m2).eq.m2))then
        //     ! upper right, m<0 + m2 > 0
        //     element = -dcmplx(0.d0,(-1.d0)**(l-(l-m)))
        // endif
        // and matmul with its T matrix
    }
    return cplxdb(0., 0.);
}

inline cplxdb get_C_matrix_element(int m, int mp, const CODE_CHOICE &choice)
{
    if (choice == CODE_CHOICE::ORIG)
        return get_C_matrix_element_orig(m, mp);
    if (choice == CODE_CHOICE::AIMS)
        return get_C_matrix_element_aims(m, mp);
    else
        throw std::invalid_argument("Unknown choice of real spherical harmonics");
}

matrix<cplxdb> get_C_matrix(unsigned l, const CODE_CHOICE &choice);
