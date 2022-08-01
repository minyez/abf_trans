#include "../src/kmesh.h"
#include "testutils.h"
#include <iostream>

using namespace std;

void test_get_kgrids()
{
    matrix<double> kgrids;

    cout << "kgrids for (2, 2, 2)" << endl;
    kgrids = get_kgrids({2, 2, 2});
    assert(kgrids.nr == 8);
    // cout << kgrids;

    cout << "kgrids for (3, 3, 3)" << endl;
    kgrids = get_kgrids({3, 3, 3});
    assert(kgrids.nr == 27);
    // cout << kgrids;

    cout << "kgrids for (4, 2, 3)" << endl;
    kgrids = get_kgrids({4, 2, 3});
    assert(kgrids.nr == 24);
    // cout << kgrids;

    cout << "kgrids for (5, 5, 5)" << endl;
    kgrids = get_kgrids({5, 5, 5});
    assert(kgrids.nr == 125);
    // cout << kgrids;
}

void test_KGrids_irkgrids_diamond()
{
    // diamond 2x2x2
    matrix<double> latt(3, 3), posi(2, 3);
    latt(0, 1) = latt(0, 2) = latt(1, 0) = latt(1, 2) = latt(2, 0) = latt(2, 1) = 0.5;
    posi(1, 0) = posi(1, 1) = posi(1, 2) = 0.25;
    SpgDS_c dataset(latt, posi, {6, 6});
    KGrids kgrids(2, 2, 2);
    kgrids.generate_irk_map(dataset);
    assert(kgrids.nirkpts == 3);
    cout << "Irreducible points for 2x2x2" << endl;
    for (int i = 0; i < kgrids.nirkpts; i++)
        cout << kgrids.irkpts.get_row(i) << " " << kgrids.irk_weights[i] << endl;

    // kgrids.rebuild_grids(3, 3, 3);
    // kgrids.generate_irk_map(dataset);
    // assert(kgrids.nirkpts == 4);
    // cout << "Irreducible points for 3x3x3" << endl;
    // cout << kgrids.irkpts;

    kgrids.rebuild_grids(4, 2, 3);
    kgrids.generate_irk_map(dataset);
    assert(kgrids.nirkpts == 12);
    cout << "Irreducible points for 4x2x3" << endl;
    for (int i = 0; i < kgrids.nirkpts; i++)
        cout << kgrids.irkpts.get_row(i) << " " << kgrids.irk_weights[i] << endl;
}

int main (int argc, char *argv[])
{
    test_get_kgrids();
    test_KGrids_irkgrids_diamond();
    return 0;
}
