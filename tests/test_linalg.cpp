#include "../src/linalg.h"
#include "testutils.h"
#include <cassert>

void test_blas_mv_mm_double(bool debug = false)
{
    const int RDIM = 6;
    const int IDIM = 2;
    const int JDIM = 3;
    assert( IDIM * JDIM == RDIM );
    const int CDIM = 2;
    double v1[RDIM] { -1, 2, 3, 1, 4, 6};
    double v2[RDIM] { 4, 0, 1, -1, 2, 3};
    double mat[RDIM][CDIM] {
        {1, -1}, {-1, 2}, {3, 0}, {0, -1}, {2, 3}, {0, 1}
    };
    int incr = 1;
    // vector dot
    assert(ddot_(&RDIM, v1, &incr, v2, &incr) == 24);
    assert(linalg::dot(RDIM, v1, incr, v2, incr) == 24);

    // matrix-vector
    double y[CDIM] {};
    const char transn = 'N';
    const double alpha = 1.0, beta = 0.0;
    dgemv_(&transn, &CDIM, &RDIM, &alpha, *mat, &CDIM, v1, &incr, &beta, y, &incr);
    assert( y[0] == 14.0 && y[1] == 22.0 );
    /* linalg::gemv('T', RDIM, CDIM, 1.0, *mat, CDIM, v1, 1, 0.5, y, 1); */
    /* assert( y[0] == 21.0 && y[1] == 33.0 ); */

    srand(time(0));

    // matrix-matrix
    double mat2[CDIM][RDIM] = {};
    for ( int ic = 0; ic != CDIM; ic++ )
        for ( int ir = 0; ir != RDIM; ir++ )
            mat2[ic][ir] = double(rand()) / RAND_MAX;

    double mat_or[RDIM][RDIM] = {};
    double mat_oc[CDIM][CDIM] = {};
    const char transt = 'T';
    // mat[RDIM][CDIM] * mat2[CDIM][RIM] = m[RDIM][RDIM]
    dgemm_(&transn, &transn, &RDIM, &RDIM, &CDIM, &alpha, *mat2, &RDIM, *mat, &CDIM, &beta, *mat_or, &RDIM);
    assert( is_matmul_AB_equal_C(RDIM, RDIM, CDIM, *mat, *mat2, *mat_or, false, debug));
    linalg::gemm(transn, transn, RDIM, RDIM, CDIM, 1.0, *mat, CDIM, *mat2, RDIM, 0.0, *mat_or, RDIM);
    assert( is_matmul_AB_equal_C(RDIM, RDIM, CDIM, *mat, *mat2, *mat_or, false, debug));
}

int main (int argc, char *argv[])
{
    test_blas_mv_mm_double(true);
    return 0;
}
