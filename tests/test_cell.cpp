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

int main (int argc, char *argv[])
{
    test_get_recip_latt();
    return 0;
}


