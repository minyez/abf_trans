#include "../src/vec.h"
#include "testutils.h"

void test_norm()
{
    double dva[] = {1.0, 2.0, 3.0};
    vec<double> dv(3, dva);
    assert(fequal(norm2(dv), sqrt(14)));
    std::complex<double> cva[] = {{1.0, 1.0}, {2.0, -1}, {3.0, -2.0}};
    vec<std::complex<double>> cv(3, cva);
    assert(fequal(norm2(cv), sqrt(20)));

    std::complex<double> cv2a[] = {{1.0, 4.0}, {2.0, -1}, {3.0, -2.0}};
    vec<std::complex<double>> cv2(3, cv2a);
    assert(cv < cv2);
}

int main (int argc, char *argv[])
{
    test_norm();
    return 0;
}
