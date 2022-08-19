#include "../src/cell.h"
#include "../src/constants.h"
#include "testutils.h"
#include <iostream>

using namespace std;

void test_get_recip_latt()
{
    matrix<double> latt(3, 3), recplatt_ref(3, 3);

    // primitive FCC as direct, BCC as reciprocal
    latt(0, 1) = latt(0, 2) = latt(1, 0) = latt(1, 2) = latt(2, 0) = latt(2, 1) = 0.5;
    latt *= 2*PI;
    recplatt_ref = 1.0;
    recplatt_ref.set_diag(-1);
    /* cout << get_recip_latt(latt) << recplatt_ref; */
    assert(get_recip_latt(latt) == recplatt_ref);

    // hexagonal
    latt.zero_out();
    recplatt_ref.zero_out();
    latt(0, 0) = 2.0;
    latt(1, 0) = -1.0;
    latt(1, 1) = std::sqrt(3);
    latt(2, 2) = 4;
    latt *= 2*PI;
    recplatt_ref(0, 0) = 0.5;
    recplatt_ref(0, 1) = 0.5 / std::sqrt(3);
    recplatt_ref(1, 1) = 1 / std::sqrt(3);
    recplatt_ref(2, 2) = 0.25;
    cout << latt << get_recip_latt(latt) << recplatt_ref;
    assert(get_recip_latt(latt) == recplatt_ref);
}

void test_move_to_center()
{
    vec<double> a(3), b(3), ref(3), shift(3);
    a[0] = 2.0, a[1] = 0.2, a[2] = -1.3;
    b = move_to_center(a, 0.0, true);
    ref[0] = 0.0, ref[1] = 0.2, ref[2] = 0.7;
    shift[0] = 2.0, shift[1] = 0.0, shift[2] = -2.0;
    assert(ref == a && shift == b);
}

int main (int argc, char *argv[])
{
    test_get_recip_latt();
    test_move_to_center();
    return 0;
}


