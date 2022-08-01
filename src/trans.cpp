#include "trans.h"
#include <set>

matrix<cplxdb> compute_W_matrix(const matrix<double> &lattice, 
                                const matrix<double> &positions, 
                                const vector<int> &atom_types,
                                const vec<double> &kprime,
                                const matrix<double> &rotmat_spg,
                                const vec<double> &transi_spg,
                                const std::map<int, std::vector<abf_id>> &map_type_abfs_in)
{
    ABF abasis(atom_types, map_type_abfs_in);
    const int nabf = abasis.get_number_of_total_abfs();
    const auto AAT = lattice * transpose(lattice);
    matrix<cplxdb> Wmat(nabf, nabf);
    // k = V k'
    const auto k = AAT * rotmat_spg * inverse(AAT) * kprime;

    for (int iMu = 0; iMu < atom_types.size(); iMu++)
    {
        const auto s_Mu = positions.get_row(iMu);
        auto sprime_Mu = inverse(rotmat_spg) * (s_Mu + transi_spg);
        const auto OVMu = move_to_center(sprime_Mu, 0.0);
        const double ang = dot(kprime, OVMu);
        const cplxdb eprefac{std::cos(ang), -std::sin(ang)};
        // debug
        std::cout << eprefac << std::endl;
        const auto abfs_iMu = map_type_abfs[atom_types[iMu]];
        vector<int> ls_compute;
        for (const auto &abf: abfs_iMu)
            ls_compute.push_back(abf.l);
        
        for (int iMup = 0; iMup < atom_types.size(); iMup++)
        {
            const auto s_Mup = positions.get_row(iMup);
            if (s_Mup == s_Mu)
            {
                // loop over l to avoid duplicate calculations of radial functions with the same l
                for (const auto &l: std::set<int>{ls_compute.cbegin(), ls_compute.cend()})
                {
                    bool is_proper;
                    auto euler = get_Euler_from_sym_matrix_spg(rotmat_spg, lattice, is_proper);
                    auto Delta = get_RSH_Delta_matrix_from_Euler(l, euler, is_proper);
                    Delta.conj();
                    std::cout << Delta; // debug
                    Delta *= eprefac;
                    for (int irf = 0; irf < abfs_iMu.size(); irf++)
                    {
                        // distribute to the full matrix
                        if (l == abfs_iMu[irf].l)
                        {
                            for (int m = -l; m <= l; m++)
                            {
                                int iabf1 = abasis.get_abf_index(iMu, irf, m);
                                for (int mp = -l; mp <= l; mp++)
                                {
                                    int iabf2 = abasis.get_abf_index(iMup, irf, mp);
                                    Wmat(iabf1, iabf2) = Delta(m+l, mp+l);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return Wmat;
}
