#include "../src/mathtools.h"
#include "testutils.h"

void test_factorials()
{
    assert(fequal(factorial(0), 1.));
    assert(fequal(factorial(1), 1.));
    assert(fequal(factorial(5), 120.));
    assert(factorial(-2) == INFINITY);

    assert(fequal(double_factorial(2), 2.));
    assert(fequal(double_factorial(-5), 1./3));
    assert(fequal(double_factorial(-9), 1./105));
    assert(double_factorial(-2) == INFINITY);
}

int main (int argc, char *argv[])
{
    test_factorials();
    return 0;
}
