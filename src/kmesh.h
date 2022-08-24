#pragma once
#include "mathtools.h"
#include "matrix.h"
#include "spglib_utils.h"
#include <array>
#include <map>

// the kgrids are always Gamma-centered
matrix<double> get_kgrids(std::array<int, 3> nks, const CODE_CHOICE &code = CODE_CHOICE::ORIG);

inline bool is_same_k(const vec<double> &k1, const vec<double> &k2, const double thres = 1.0e-5)
{
    auto kdiff = k1 - k2;
    for (int i = 0; i < kdiff.size(); i++)
    {
        shift_to_unit(kdiff[i], 0.0, true);
        if(1 - kdiff[i] < thres) kdiff[i] = 1 - kdiff[i];
    }
    return vec_equal(kdiff, {3}, thres);
}

class KGrids
{
private:
    void set_kgrids();
    void clear_irkgrids();
public:
    CODE_CHOICE code;
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
    //! like indexmap_k_symop, but store the symmetry operation of inverse mapping, i.e. from irk to k
    std::vector<int> indexmap_k_invsymop;
    
    std::map<int, std::vector<int>> map_iirk_iks;

    // construtors and destructors
    KGrids(int nkx, int nky, int nkz, const CODE_CHOICE &code_in = CODE_CHOICE::ORIG);
    KGrids(const std::array<int, 3> &nks_in, const CODE_CHOICE &code_in = CODE_CHOICE::ORIG);
    ~KGrids() {};
    void rebuild_grids(int nkx, int nky, int nkz, const CODE_CHOICE &code_in);
    void rebuild_grids(const std::array<int, 3> &nks_in, const CODE_CHOICE &code_in);
    int index(const vec<double> &kvec) const
    {
        // a naive implementation of find
        for (int i = 0; i < kpts.nr; i++)
            if (is_same_k(kpts.get_row(i), kvec, 1e-6)) return i;
        return -1;
    }
    int index(double k1, double k2, double k3) const
    {
        double k[3] {k1, k2, k3};
        return index({3, k});
    }
    int index(const double *kvec) const { return index({3, kvec}); }
    bool have_k(const vec<double> &kvec) const { return index(kvec) != -1; }
    bool have_k(const double *kvec) const { return index(kvec) != -1; }
    bool have_k(double k1, double k2, double k3) const { return index(k1, k2, k3) != -1; }

    int index_irk(const vec<double> &kvec) const
    {
        if (!irkgrids_generated()) throw std::logic_error("irkgrids not generated");
        // a naive implementation of find
        for (int i = 0; i < irkpts.nr; i++)
            if (is_same_k(irkpts.get_row(i), kvec, 1e-6)) return i;
        return -1;
    }
    bool have_irk(const vec<double> &kvec) const { return index_irk(kvec) != -1; }

    // member functions
    void generate_irk_map(const SpgDS_c &dataset);
    bool irkgrids_generated() const { return !indexmap_k_irk.empty(); }
};

// k in reciprocal lattice vector coordinate
void get_all_equiv_k(const vec<double> &k,
                     const vector<matrix<int>> &rotmats_spg,
                     vector<vec<double>> &equiv_ks, vector<int> &irots,
                     const CODE_CHOICE &code);

//! get all symmetry operations (index) V in rotmats_spg such that Vk1 = k2
vector<int> get_symops_connecting_k1_k2(const vec<double> &k1, const vec<double> &k2,
                                        const vector<matrix<int>> &rotmats_spg);

matrix<double> move_k_back(matrix<double> &kpts, const CODE_CHOICE &code);
vec<double> move_k_back(vec<double> &kpts, const CODE_CHOICE &code);
