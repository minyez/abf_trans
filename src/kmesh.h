#pragma once
#include "matrix.h"
#include "spglib_utils.h"
#include <map>
#include <utility>
#include <array>

// the kgrids are always Gamma-centered
matrix<double> get_kgrids(std::array<int, 3> nks);

class KGrids
{
public:
    // member attributes
    std::array<int, 3> nks;
    matrix<double> kpts;
    matrix<double> irkpts;
    std::vector<std::pair<int, int>> pair_k_irk;
    std::map<std::pair<int, int>, int> map_kirkpair_symop;

    // construtors and destructors
    KGrids(int nkx, int nky, int nkz);
    KGrids(std::array<int, 3> nks_in);
    ~KGrids() {};

    // member functions
    void generate_irk_map(const SpgDS_c &dataset);
    bool irk_map_generated() { return !pair_k_irk.empty(); }
};
