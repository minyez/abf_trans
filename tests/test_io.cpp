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
        cout << "decoding fraction string " << frac_str << ", get " << a << endl;
        assert(fequal(a, fracs[i]));
    }
}

int main (int argc, char *argv[])
{
    test_decode_fraction();
    return 0;
}
