#include <iostream>
#include <cassert>
#include "../src/sph.h"
#include "testutils.h"

void test_C_matrix_l1()
{
    const int l = 1;
    const unsigned msize = 2*l + 1;

    auto cmat = get_C_matrix(l, CODE_CHOICE::ORIG);
    matrix<cplxdb> cmat_ref(msize, msize);

    // initialize
    cmat_ref(0, 0) = cmat_ref(0, 2) = cplxdb{0, 1};
    cmat_ref(1, 1) = cplxdb{std::sqrt(2), 0};
    cmat_ref(2, 0) = cplxdb{1, 0};
    cmat_ref(2, 2) = cplxdb{-1, 0};

    cmat_ref /= std::sqrt(2);

    assert(is_mat_A_equal_B(msize, msize, cmat.c, cmat_ref.c, false, true));
}


int main (int argc, char *argv[])
{
    test_C_matrix_l1();
    return 0;
}
