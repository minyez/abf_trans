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
private:
    void set_kgrids();
    void clear_irkgrids();
public:
    // member attributes
    std::array<int, 3> nks;
    //! all k-points
    matrix<double> kpts;
    int nkpts;
    //! irreducible k-points
    matrix<double> irkpts;
    int nirkpts;
    //! index of irreducible k in all k-points (kpts), size = irkpts.size()
    std::vector<int> irk_index;
    //! integer weights of irreducible k, size = irkpts.size()
    std::vector<int> irk_weights;
    //! map of k-point index to its equivalent irreducible kpoints, size = kpts.size()
    std::vector<int> indexmap_k_irk;
    //! map of k-point index to the symmetry operation to get its equivalent irreducible kpoints, size = kpts.size()
    std::vector<int> indexmap_k_symop;

    // construtors and destructors
    KGrids(int nkx, int nky, int nkz);
    KGrids(const std::array<int, 3> &nks_in);
    ~KGrids() {};
    void rebuild_grids(int nkx, int nky, int nkz);
    void rebuild_grids(const std::array<int, 3> &nks_in);

    // member functions
    void generate_irk_map(const SpgDS_c &dataset);
    bool irkgrids_generated() { return !indexmap_k_irk.empty(); }
};
