#include "constants.h"
#include "cell.h"

matrix<double> latt;
matrix<double> posi_frac;
vector<int> types;
set<int> inequiv_types;
map<int, vector<int>> map_type_iatoms;

// both a and b vectors are aligned by rows
matrix<double> get_recip_latt(const matrix<double> &latt)
{
    auto recip_latt = inverse(latt) * 2.0 * PI;
    return transpose(recip_latt);
}

void generate_map_type_iatom(const vector<int> &atypes, set<int> &inequiv_types, map<int, vector<int>> &amap)
{
    amap.clear();
    inequiv_types.clear();

    inequiv_types = set<int>(atypes.cbegin(), atypes.cend());

    for (int i = 0; i < atypes.size(); i++)
        amap[atypes[i]].push_back(i);
}

matrix<double> move_to_center_cell(matrix<double> &posi, const double lowlim)
{
    matrix<double> R = posi;
    R.zero_out();
    const double uplim = lowlim + 1.0;
    for (int ia = 0; ia < R.nr; ia++)
        for (int ic = 0; ic < R.nc; ic++)
        {
            while (posi(ia, ic) < lowlim)
            {
                posi(ia, ic) += 1.0;
                R(ia, ic) -= 1.0;
            }
            while (posi(ia, ic) >= uplim)
            {
                posi(ia, ic) -= 1.0;
                R(ia, ic) += 1.0;
            }
        }
    return R;
}

vec<double> move_to_center_cell(vec<double> &posi, const double lowlim)
{
    vec<double> R = posi;
    const double uplim = lowlim + 1.0;
    R.zero_out();
    for (int ic = 0; ic < R.n; ic++)
    {
        while (posi[ic] < lowlim)
        {
            posi[ic] += 1.0;
            R[ic] -= 1.0;
        }
        while (posi[ic] >= uplim)
        {
            posi[ic] -= 1.0;
            R[ic] += 1.0;
        }
    }
    return R;
}
