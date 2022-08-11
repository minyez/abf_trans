#include "../src/io.h"
#include "testutils.h"
#include <iostream>

using namespace std;

void test_decode_fraction()
{
    double a;
    const vector<string> frac_strs = {"0.5", "1/2", "1/3", "5 / 8", "-2.0 / 0.5"};
    const vector<double> fracs = {0.5, 0.5, 1.0/3.0, 5.0/8.0, -4.0};

    assert(frac_strs.size() == fracs.size());

    for (int i = 0; i < frac_strs.size(); i++)
    {
        string frac_str = frac_strs[i];
        a = decode_fraction(frac_str);
        // cout << "decoding fraction string " << frac_str << ", get " << a << endl;
        assert(fequal(a, fracs[i]));
    }
    printf("test of decoding fraction passed\n");
}

void test_read_csc_complex_matrix()
{
    const auto cmat = read_csc_cplxdb("sample_complexmat.csc");
    const vector<int> rows
        {1, 136, 30, 22, 42, 128, 23, 149, 29, 63, 44, 258, 230};
    const vector<int> cols
        {1, 2, 6, 8, 10, 12, 14, 15, 20, 33, 50, 62, 65};
    const vector<cplxdb> vals
        {
            { 9.999823902394842e-01, 0.000000000000000e+00},
            { 1.781588226665285e-04,-1.028600442250274e-04},
            { 1.588186776101813e-22, 3.310976225044609e-06},
            { 0.000000000000000e+00,-2.165905580137130e-07},
            { 3.899897152846147e-05, 7.321788384219558e-22},
            { 2.165824469666926e-08, 3.751318021750593e-08},
            { 8.077935669463161e-28, 1.237317758123938e-10},
            {-2.269780979640086e-07,-3.931375978817807e-07},
            { 2.944303070553680e-05, 1.864264440921041e-22},
            { 2.201507879500560e-02, 5.886878983417387e-20},
            { 1.052152243175773e-02,-3.498236800428520e-19},
            { 1.366108858851323e-02, 2.366169952199770e-02},
            { 5.212060608805230e-02,-3.009184595526304e-02},
        };
    assert(rows.size() == cols.size() && rows.size() == vals.size());
    for (int i = 0; i < rows.size(); i++)
    {
        assert(fequal(vals[i], cmat(rows[i]-1, cols[i]-1), 1e-10));
    }
    printf("test of reading CSC complex matrix passed\n");

// 128 12 2.165824469666926e-08 3.751318021750593e-08
// 23 14 8.077935669463161e-28 1.237317758123938e-10
// 149 15 -2.269780979640086e-07 -3.931375978817807e-07
// 29 20 2.944303070553680e-05 1.864264440921041e-22
// 63 33 2.201507879500560e-02 5.886878983417387e-20
// 44 50 1.052152243175773e-02 -3.498236800428520e-19
// 258 62 1.366108858851323e-02 2.366169952199770e-02
// 230 65 5.212060608805230e-02 -3.009184595526304e-02
}

int main (int argc, char *argv[])
{
    test_decode_fraction();
    test_read_csc_complex_matrix();
    return 0;
}
