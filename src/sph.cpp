#include <cmath>
#include "sph.h"

matrix<cplxdb> get_C_matrix(unsigned l)
{
    auto msize = get_msize(l);
    const int il = int(l);
    matrix<cplxdb> cmat(msize, msize);
    for (int m = -il; m <= il; m++)
        for (int mp = -il; mp <= il; mp++)
        {
            cmat(m+il, mp+il) = get_C_matrix_element(m, mp);
            /* std::cout << m << " " << mp << " " << get_C_matrix_element(l, m, mp) << std::endl; */
        }
    return cmat;
}

matrix<cplxdb> get_C_matrix_aims(unsigned l)
{
    auto msize = get_msize(l);
    const int il = int(l);
    matrix<cplxdb> cmat(msize, msize);
    for (int m = -il; m <= il; m++)
        for (int mp = -il; mp <= il; mp++)
        {
            cmat(m+il, mp+il) = get_C_matrix_element_aims(m, mp);
            /* std::cout << m << " " << mp << " " << get_C_matrix_element(l, m, mp) << std::endl; */
        }
    return cmat;
}

