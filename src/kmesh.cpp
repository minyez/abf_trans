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

    return matrix<double>(kgrids_unsorted);
}

KGrids::KGrids(int nkx, int nky, int nkz): nks({nkx, nky, nkz})
{
    kpts = get_kgrids(nks);
}

KGrids::KGrids(std::array<int, 3> nks_in): nks(nks_in)
{
    kpts = get_kgrids(nks);
}

void KGrids::generate_irk_map(const SpgDS_c &dataset)
{
    vector<vec<double>> irkpts_vec;
    irkpts = irkpts_vec;
    for (int ik = 0; ik < kpts.nr; ik++)
    {
        const auto & kpt = kpts.get_row(ik);
    }
}
