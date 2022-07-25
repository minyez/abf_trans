#include <iostream>
#include <cassert>
#include "../src/sph.h"
#include "testutils.h"

void test_C_matrix_l1()
{
    const int l = 1;
    const unsigned msize = 2*l + 1;

    cplxdb *cmat = new cplxdb [msize*msize];
    get_C_matrix(cmat, l);

    cplxdb *cmat_ref = new cplxdb [msize*msize];

    // initialize
    for (int i = 0; i != msize*msize; i++)
        cmat_ref[i] = 0;
    cmat_ref[0] = cmat_ref[2] = cplxdb(0, 1);
    cmat_ref[4] = cplxdb(std::sqrt(2), 0);
    cmat_ref[6] = cplxdb(1, 0);
    cmat_ref[8] = cplxdb(-1, 0);

    for (int i = 0; i != msize*msize; i++)
        cmat_ref[i] /= std::sqrt(2);

    assert(is_mat_A_equal_B(msize, msize, cmat, cmat_ref, false, true));

    delete [] cmat;
    delete [] cmat_ref;
}


int main (int argc, char *argv[])
{
    test_C_matrix_l1();
    return 0;
}
