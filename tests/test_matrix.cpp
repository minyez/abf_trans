#include "../src/matrix.h"
#include "testutils.h"
#include <cassert>
#include <iostream>

void test_add_subtract()
{
    using namespace std;
    matrix<double> a(2, 2);
    a(0, 0) = 2;
    a(1, 1) = 1;
    a += 1;
    assert(fequal(a(0, 0), 3.0) && fequal(a(1, 1), 2.0));

    a(0, 1) = 0;
    auto aa = a * a;
    cout << "a * a: " << endl;
    cout << "  " << aa(0, 0) << " " << aa(0, 1) << endl;
    cout << "  " << aa(1, 0) << " " << aa(1, 1) << endl;
    aa -= 5;
    assert(fequal(aa(0, 0), 4.0) && fequal(aa(0, 1), -5.0) &&
           fequal(aa(1, 0), 0.0) && fequal(aa(1, 1), -1.0));
}

void test_inverse()
{
    using namespace std;
    matrix<double> a(2, 2);
    a(0, 0) = 0.5;
    a(0, 1) = 2;
    a(1, 1) = 1;
    auto b = inverse(a);
    cout << b.size() << endl;
    cout << "inverse of a: " << endl;
    cout << "  " << b(0, 0) << " " << b(0, 1) << endl;
    cout << "  " << b(1, 0) << " " << b(1, 1) << endl;
    assert(fequal(b(0, 0), 2.0) && fequal(b(0, 1), -4.0) &&
           fequal(b(1, 0), 0.0) && fequal(b(1, 1), 1.0));
}

int main (int argc, char *argv[])
{
    test_add_subtract();
    test_inverse();

    return 0;
}
