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
    cout << kgrids;

    cout << "kgrids for (3, 3, 3)" << endl;
    kgrids = get_kgrids({3, 3, 3});
    assert(kgrids.nr == 27);
    cout << kgrids;

    cout << "kgrids for (5, 5, 5)" << endl;
    kgrids = get_kgrids({5, 5, 5});
    assert(kgrids.nr == 125);
    cout << kgrids;
}

int main (int argc, char *argv[])
{
    test_get_kgrids();
    return 0;
}
