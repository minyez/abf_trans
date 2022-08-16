#include "trans.h"
#include "constants.h"
#include <set>
#include <iostream>

matrix<cplxdb> compute_W_matrix(const matrix<double> &lattice, 
                                const matrix<double> &positions, 
                                const vector<int> &atom_types,
                                const vec<double> &kprime,
                                const matrix<int> &rotmat_spg,
                                const vec<double> &transi_spg,
                                const std::map<int, std::vector<abf_id>> &map_type_abfs_in,
                                const CODE_CHOICE & choice)
{
    ABF abasis(atom_types, map_type_abfs_in);
    const int nabf = abasis.get_number_of_total_abfs();
    const auto AAT = lattice * transpose(lattice);
    matrix<cplxdb> Wmat(nabf, nabf);
    // const double prefac = choice == CODE_CHOICE::AIMS ? -1. : 1.;
    const double prefac = choice == CODE_CHOICE::AIMS ? -1. : 1.;
    const auto rotmat_spg_db = to_double(rotmat_spg);
    // k = V k'
    const auto k = AAT * rotmat_spg_db * inverse(AAT) * kprime;

    for (int iMu = 0; iMu < atom_types.size(); iMu++)
    {
        const auto s_Mu = positions.get_row(iMu);
        auto sprime_Mu = inverse(rotmat_spg_db) * (s_Mu - transi_spg);
        // std::cout << "S'_Mu = " << sprime_Mu << std::endl; // debug
        const auto OVMu = move_to_center(sprime_Mu, 0.0, true);
        // std::cout << "Shift back, S'_Mu = " << sprime_Mu << " , OV_Mu = " << OVMu << std::endl; // debug
        const double ang = prefac * 2.0 * PI * dot(kprime, OVMu);
        const cplxdb eprefac(std::cos(ang), -std::sin(ang));
        // std::cout << ang << " " << eprefac << std::endl; // debug
        const auto abfs_iMu = map_type_abfs[atom_types[iMu]];
        vector<int> ls_compute;
        for (const auto &abf: abfs_iMu)
            ls_compute.push_back(abf.l);

        for (int iMup = 0; iMup < atom_types.size(); iMup++)
        {
            const auto s_Mup = positions.get_row(iMup);
            if (s_Mup == sprime_Mu)
            {
                // loop over l to avoid duplicate calculations of radial functions with the same l
                for (const auto &l: std::set<int>{ls_compute.cbegin(), ls_compute.cend()})
                {
                    bool is_proper;
                    auto euler = get_Euler_from_sym_matrix_spg(rotmat_spg, lattice, is_proper);
                    auto Delta = get_RSH_Delta_matrix_from_Euler(l, euler, is_proper, choice);
                    Delta.conj();
                    // std::cout << Delta; // debug
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
                                    // printf("%5d %5d %5d %5d %5d %5d %5d %5d %10.5f %10.5f\n", iMu, iMup, irf, l, m, mp, iabf1, iabf2, Wmat(iabf1, iabf2).real(), Wmat(iabf1, iabf2).imag()); // debug
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

matrix<cplxdb> compute_representation_on_equiv_k(const vec<double> &kprime, const matrix<cplxdb> &cmat_at_kprime,
                                                 const matrix<double> &lattice,
                                                 const matrix<double> &positions,
                                                 const vector<int> &atom_types,
                                                 const matrix<int> &rotmat_spg, const vec<double> &transi_spg,
                                                 const std::map<int, std::vector<abf_id>> &map_type_abfs_in,
                                                 const CODE_CHOICE & choice)
{
    auto wmat = compute_W_matrix(lattice, positions, atom_types, kprime, rotmat_spg, transi_spg, map_type_abfs_in, choice);
    return wmat * cmat_at_kprime * transpose(wmat, true);
}
