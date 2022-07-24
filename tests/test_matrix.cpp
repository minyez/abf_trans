#include "../src/matrix.h"
#include "testutils.h"
#include <cassert>
#include <iostream>
#include <vector>

void test_add_subtract()
{
    using namespace std;
    matrix<double> a(2, 2);
    a(0, 0) = 2;
    a(1, 1) = 1;
    a += 1;
    assert(fequal(a(0, 0), 3.0) && fequal(a(1, 1), 2.0));

    a(0, 1) = 0;
    matrix<double> aa = a * a;
    aa -= 5;
    cout << "a * a - 5: " << endl;
    cout << "  " << aa(0, 0) << " " << aa(0, 1) << endl;
    cout << "  " << aa(1, 0) << " " << aa(1, 1) << endl;
    assert(fequal(aa(0, 0), 4.0) && fequal(aa(0, 1), -5.0) &&
           fequal(aa(1, 0), 0.0) && fequal(aa(1, 1), -1.0));

    double arr[6] = {1, 2, 3, 4, 5, 6};
    matrix<double> b(3, 2, arr);
    vector<double> v = {1, -1};
    auto bplusv = b + v;
    cout << "b: " << endl;
    cout << "  " << b(0, 0) << " " << b(0, 1) << endl;
    cout << "  " << b(1, 0) << " " << b(1, 1) << endl;
    cout << "  " << b(2, 0) << " " << b(2, 1) << endl;
    cout << "b + v: " << endl;
    cout << "  " << bplusv(0, 0) << " " << bplusv(0, 1) << endl;
    cout << "  " << bplusv(1, 0) << " " << bplusv(1, 1) << endl;
    cout << "  " << bplusv(2, 0) << " " << bplusv(2, 1) << endl;
    assert(fequal(bplusv(0, 0), 2.0) && fequal(bplusv(0, 1), 1.0) &&
           fequal(bplusv(1, 0), 4.0) && fequal(bplusv(1, 1), 3.0) &&
           fequal(bplusv(2, 0), 6.0) && fequal(bplusv(2, 1), 5.0));
}

void test_inverse_det()
{
    using namespace std;
    matrix<double> a(2, 2);
    a(0, 0) = 0.5;
    a(0, 1) = 2;
    a(1, 1) = 1;
    auto b = inverse(a);
    cout << b.size() << endl;
    cout << "inverse of a: " << endl;
    cout << b;
    assert(fequal(b(0, 0), 2.0) && fequal(b(0, 1), -4.0) &&
           fequal(b(1, 0), 0.0) && fequal(b(1, 1), 1.0));
    cout << "Determinant = " << b.det() << endl;
    assert(b.det() == 2);
    b += 1;
    cout << "add 1 to a^-1: " << endl;
    cout << b;
    cout << "Determinant = " << b.det() << endl;
    assert(b.det() == 9);
    b -= 4;
    cout << "then subtract 4" << endl;
    cout << b;
    cout << "Determinant = " << b.det() << endl;
    assert(b.det() == -19);
}

int main(int argc, char *argv[])
{
    test_add_subtract();
    test_inverse_det();

    return 0;
}
