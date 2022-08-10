#pragma once
#include "matrix.h"
#include "spglib_utils.h"
#include <array>

// the kgrids are always Gamma-centered
matrix<double> get_kgrids(std::array<int, 3> nks, const CODE_CHOICE &code = CODE_CHOICE::ORIG);

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
            if (vec_equal(kpts.get_row(i), kvec, 1e-6)) return i;
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

    // member functions
    void generate_irk_map(const SpgDS_c &dataset);
    bool irkgrids_generated() { return !indexmap_k_irk.empty(); }
};

// k in reciprocal lattice vector coordinate
void get_all_equiv_k(const vec<double> &k, const matrix<double> lattice,
                     const vector<matrix<int>> &rotmats_spg,
                     vector<vec<double>> &equiv_ks, vector<int> &irots,
                     const CODE_CHOICE &code);

//! get all symmetry operations (index) V in rotmats_spg connecting k to Vk
vector<int> get_all_symops_connecting_ks(const vec<double> &k, const vec<double> &Vk,
                                         const matrix<double> lattice,
                                         const vector<matrix<int>> &rotmats_spg);

bool is_same_k(const vec<double> &k1, const vec<double> &k2, const double thres = 1.0e-5);

matrix<double> move_k_back(matrix<double> &kpts, const CODE_CHOICE &code);
vec<double> move_k_back(vec<double> &kpts, const CODE_CHOICE &code);
