#include <cmath>
#include "sph.h"

void get_C_matrix(cplxdb *cmat, unsigned l)
{
    auto msize = get_msize(l);
    const int il = l;
    for (int m = -l; m <= il; m++)
        for (int mp = -l; mp <= il; mp++)
        {
            cmat[(l+m)*msize+l+mp] = get_C_matrix_element(l, m, mp);
            /* std::cout << m << " " << mp << " " << get_C_matrix_element(l, m, mp) << std::endl; */
        }
}
