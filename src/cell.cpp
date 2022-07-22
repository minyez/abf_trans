#include "constants.h"
#include "cell.h"

matrix<double> latt;
matrix<double> posi_frac;
vector<int> types;
set<int> inequiv_types;
map<int, vector<int>> map_type_iatoms;

inline matrix<double> get_recip_latt(matrix<double> latt)
{
    return inverse(latt) * 2.0 * PI;
}

void generate_map_type_iatom(const vector<int> &atypes, set<int> &inequiv_types, map<int, vector<int>> &amap)
{
    amap.clear();
    inequiv_types.clear();

    inequiv_types = set<int>(atypes.cbegin(), atypes.cend());

    for (int i = 0; i < atypes.size(); i++)
        amap[atypes[i]].push_back(i);
}
