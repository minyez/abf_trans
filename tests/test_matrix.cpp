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

void test_mm_mv()
{
    using namespace std;
    double arr_a[] = {1, -1, 2, 0};
    matrix<double> a(2, 2, arr_a);
    vec<double> v(2);
    v[0] = 1;
    v[1] = -1;
    vec<double> ref(2);

    ref[0] = ref[1] = 2;
    assert(a * v == ref);

    ref[0] = ref[1] = -1;
    assert(v * a == ref);

    cplxdb arr_ca[] = {{-1.0, 0.0}, {0.0, 1.0}, {0.0, 2.0}, {2.0, 0.0}};
    matrix<cplxdb> ca(2, 2, arr_ca);
    cplxdb arr_cv[] = {{1.0, 0.0}, {0.0, 1.0}};
    vec<cplxdb> cv(2, arr_cv);

    vec<cplxdb> cref(2);
    cref[0] = cplxdb{-2.0, 0.0};
    cref[1] = cplxdb{0.0, 4.0};
    cout << "Computing complex matrix x complex vec" << endl;
    cout << ca * cv;
    assert(ca * cv == cref);
    cref[0] = cplxdb{-3.0, 0.0};
    cref[1] = cplxdb{0.0, 3.0};
    cout << "Computing complex vec x complex matrix" << endl;
    cout << cv * ca;
    assert(cv * ca == cref);
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

    a.resize(3, 3);
    a(0, 0) = a(0, 1) = a(1, 1) = a(2, 2) = 1;
    matrix<double> aref(3, 3);
    aref(0, 0) = aref(1, 1) = aref(2, 2) = 1;
    aref(0, 1) = -1;
    cout << inverse(a) << aref;
    assert(inverse(a) == aref);
}

int main(int argc, char *argv[])
{
    test_add_subtract();
    test_mm_mv();
    test_inverse_det();

    return 0;
}
