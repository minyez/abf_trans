#include "kmesh.h"
#include <algorithm>
#include <vector>

matrix<double> get_kgrids(std::array<int, 3> nks)
{
    const int nkpts = nks[0] * nks[1] * nks[2];
    std::vector<vec<double>> kgrids_unsorted(nkpts);
    int ikpt = 0;
    for (int ikx = 0; ikx < nks[0]; ikx++)
    {
        double kx = double(ikx) / nks[0];
        for (int iky = 0; iky < nks[1]; iky++)
        {
            double ky = double(iky) / nks[1];
            for (int ikz = 0; ikz < nks[2]; ikz++)
            {
                double kz = double(ikz) / nks[2];
                kgrids_unsorted[ikpt].resize(3);
                kgrids_unsorted[ikpt][0] = kx;
                kgrids_unsorted[ikpt][1] = ky;
                kgrids_unsorted[ikpt][2] = kz;
                ikpt++;
            }
        }
    }
    vec<double> shift(3);
    for (int i = 0; i < 3; i++)
    {
        shift[0] = ( nks[i]%2 ? nks[i] / 2 : (nks[i] - 1) / 2 ) / double(nks[i]);
    }
    for (int i = 0; i < nkpts; i++)
        kgrids_unsorted[i] -= shift;

    sort(kgrids_unsorted.begin(), kgrids_unsorted.end());
    matrix<double> kgrids(nkpts, 3);
    for (int ikpt = 0; ikpt < nkpts; ikpt++)
        for (int i = 0; i < 3; i++)
            kgrids(ikpt, i) = kgrids_unsorted[ikpt][i];

    return kgrids;
}
