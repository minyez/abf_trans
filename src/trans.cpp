#include "trans.h"
#include "constants.h"
#include "kmesh.h"
#include <set>
#include <iostream>

matrix<cplxdb> compute_M_matrix(const matrix<double> &lattice, 
                                const matrix<double> &positions, 
                                const vector<int> &atom_types,
                                const vec<double> &k,
                                const matrix<int> &rotmat_spg,
                                const vec<double> &transi_spg,
                                const std::map<int, std::vector<abf_id>> &map_type_abfs_in,
                                const CODE_CHOICE & choice)
{
    ABF basis(atom_types, map_type_abfs_in);
    const int nabf = basis.get_number_of_total_abfs();
    int a, b;
    get_bloch_phase_convention(choice, a, b);
    matrix<cplxdb> Mmat(nabf, nabf);
    const auto rotmat_spg_db = to_double(rotmat_spg);
    // ~k = V k
    auto tilde_k = inverse(transpose(rotmat_spg_db)) * k;
    move_k_back(tilde_k, choice);

    auto rotmat_xyz = get_sym_matrix_xyz(rotmat_spg, lattice);
    bool is_proper;
    auto euler = get_Euler_from_sym_matrix_xyz(inverse(rotmat_xyz), is_proper);

    for (int iMu = 0; iMu < atom_types.size(); iMu++)
    {
        const auto s_Mu = positions.get_row(iMu);
        const auto s_tildeMu = rotmat_spg_db * s_Mu + transi_spg;
        const auto abfs_iMu = map_type_abfs[atom_types[iMu]];
        vector<int> ls_compute;
        for (const auto &abf: abfs_iMu)
            ls_compute.push_back(abf.l);

        for (int iMup = 0; iMup < atom_types.size(); iMup++)
        {
            const auto s_Mup = positions.get_row(iMup);
            if (is_same_atom_in_center(s_Mup, s_tildeMu))
            {
                const auto OSMu = s_tildeMu - s_Mup;
                // loop over l to avoid duplicate calculations of radial functions with the same l
                for (const auto &l: std::set<int>{ls_compute.cbegin(), ls_compute.cend()})
                {
                    // auto euler = get_Euler_from_sym_matrix_spg(rotmat_spg, lattice, is_proper);
                    auto Delta = get_RSH_Delta_matrix_from_Euler(l, euler, is_proper, choice);
                    double ang = dot(tilde_k, OSMu);
                    if (b != 0)
                        ang += b * (dot(tilde_k, s_Mup) - dot(k, s_Mu));
                    ang *= a * 2.0 * PI;
                    // perform conjugate is just a formal treatment and should not affect the result, as Delta is real by definition
                    // std::cout << Delta; // debug
                    const cplxdb phase(std::cos(ang), std::sin(ang));
                    Delta *= phase;
                    for (int irf = 0; irf < abfs_iMu.size(); irf++)
                    {
                        // distribute to the full matrix
                        if (l == abfs_iMu[irf].l)
                        {
                            for (int m = -l; m <= l; m++)
                            {
                                int iabf1 = basis.get_abf_index(iMu, irf, m);
                                for (int mp = -l; mp <= l; mp++)
                                {
                                    int iabf2 = basis.get_abf_index(iMup, irf, mp);
                                    Mmat(iabf1, iabf2) = Delta(m+l, mp+l);
                                    // printf("%5d %5d %5d %5d %5d %5d %5d %5d %10.5f %10.5f\n", iMu, iMup, irf, l, m, mp, iabf1, iabf2, Wmat(iabf1, iabf2).real(), Wmat(iabf1, iabf2).imag()); // debug
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return Mmat;
}

matrix<cplxdb> compute_M_matrix_aims(const matrix<double> &lattice,
                                     const matrix<double> &positions,
                                     const vector<int> &atom_types,
                                     const vec<double> &k,
                                     const matrix<int> &rotmat_spg, const vec<double> &transi_spg,
                                     const std::map<int, std::vector<abf_id>> &map_type_abfs_in)
{
    ABF basis(atom_types, map_type_abfs_in);
    const int nabf = basis.get_number_of_total_abfs();
    matrix<cplxdb> Mmat(nabf, nabf);
    const auto rotmat_spg_db = to_double(rotmat_spg);
    auto Vk = inverse(transpose(rotmat_spg_db)) * k;
    // std::cout << k << " " << Vk << std::endl;

    auto rotmat_xyz = get_sym_matrix_xyz(rotmat_spg, lattice);
    bool is_proper;
    auto euler = get_Euler_from_sym_matrix_xyz(inverse(rotmat_xyz), is_proper);

    for (int iMu = 0; iMu < atom_types.size(); iMu++)
    {
        const auto s_Mu = positions.get_row(iMu);
        auto s_tildeMu = inverse(rotmat_spg_db) * (s_Mu - transi_spg);
        const auto abfs_iMu = map_type_abfs[atom_types[iMu]];
        vector<int> ls_compute;
        for (const auto &abf: abfs_iMu)
            ls_compute.push_back(abf.l);

        for (int iMup = 0; iMup < atom_types.size(); iMup++)
        {
            const auto s_Mup = positions.get_row(iMup);
            if (is_same_atom_in_center(s_Mup, s_tildeMu))
            {
                // loop over l to avoid duplicate calculations of radial functions with the same l
                for (const auto &l: std::set<int>{ls_compute.cbegin(), ls_compute.cend()})
                {
                    // auto euler = get_Euler_from_sym_matrix_spg(rotmat_spg, lattice, is_proper);
                    auto Delta = get_RSH_Delta_matrix_from_Euler(l, euler, is_proper, CODE_CHOICE::AIMS);
                    double ang = 2.0*PI*(dot(Vk, s_Mu) - dot(k, s_Mup));
                    const cplxdb phase(std::cos(ang), std::sin(ang));
                    Delta *= phase;
                    for (int irf = 0; irf < abfs_iMu.size(); irf++)
                    {
                        // distribute to the full matrix
                        if (l == abfs_iMu[irf].l)
                        {
                            for (int m = -l; m <= l; m++)
                            {
                                int iabf1 = basis.get_abf_index(iMu, irf, m);
                                for (int mp = -l; mp <= l; mp++)
                                {
                                    int iabf2 = basis.get_abf_index(iMup, irf, mp);
                                    Mmat(iabf1, iabf2) = Delta(m+l, mp+l);
                                    // printf("%5d %5d %5d %5d %5d %5d %5d %5d %10.5f %10.5f\n", iMu, iMup, irf, l, m, mp, iabf1, iabf2, Wmat(iabf1, iabf2).real(), Wmat(iabf1, iabf2).imag()); // debug
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return Mmat;
}

matrix<cplxdb> compute_M_matrix_abacus(const matrix<double> &lattice,
                                       const matrix<double> &positions,
                                       const vector<int> &atom_types,
                                       const vec<double> &k,
                                       const matrix<int> &rotmat_spg, const vec<double> &transi_spg,
                                       const std::map<int, std::vector<abf_id>> &map_type_abfs_in)
{
    ABF basis(atom_types, map_type_abfs_in);
    const int nabf = basis.get_number_of_total_abfs();
    matrix<cplxdb> Mmat(nabf, nabf);
    const auto rotmat_spg_db = to_double(rotmat_spg);
    auto Vk = inverse(transpose(rotmat_spg_db)) * k;
    // std::cout << k << " " << Vk << std::endl;

    auto rotmat_xyz = get_sym_matrix_xyz(rotmat_spg, lattice);
    bool is_proper;
    auto euler = get_Euler_from_sym_matrix_xyz(inverse(rotmat_xyz), is_proper);

    for (int iMu = 0; iMu < atom_types.size(); iMu++)
    {
        const auto s_Mu = positions.get_row(iMu);
        auto s_tildeMu = inverse(rotmat_spg_db) * (s_Mu - transi_spg);
        const auto abfs_iMu = map_type_abfs[atom_types[iMu]];
        vector<int> ls_compute;
        for (const auto &abf: abfs_iMu)
            ls_compute.push_back(abf.l);

        for (int iMup = 0; iMup < atom_types.size(); iMup++)
        {
            const auto s_Mup = positions.get_row(iMup);
            if (is_same_atom_in_center(s_Mup, s_tildeMu))
            {
                // loop over l to avoid duplicate calculations of radial functions with the same l
                for (const auto &l: std::set<int>{ls_compute.cbegin(), ls_compute.cend()})
                {
                    // auto euler = get_Euler_from_sym_matrix_spg(rotmat_spg, lattice, is_proper);
                    auto Delta = get_RSH_Delta_matrix_from_Euler(l, euler, is_proper, CODE_CHOICE::ABACUS);
                    double ang = -2.0*PI*(dot(Vk, s_Mu) - dot(k, s_Mup));
                    const cplxdb phase(std::cos(ang), std::sin(ang));
                    Delta *= phase;
                    for (int irf = 0; irf < abfs_iMu.size(); irf++)
                    {
                        // distribute to the full matrix
                        if (l == abfs_iMu[irf].l)
                        {
                            for (int m = -l; m <= l; m++)
                            {
                                int iabf1 = basis.get_abf_index(iMu, irf, m);
                                for (int mp = -l; mp <= l; mp++)
                                {
                                    int iabf2 = basis.get_abf_index(iMup, irf, mp);
                                    Mmat(iabf1, iabf2) = Delta(m+l, mp+l);
                                    // printf("%5d %5d %5d %5d %5d %5d %5d %5d %10.5f %10.5f\n", iMu, iMup, irf, l, m, mp, iabf1, iabf2, Wmat(iabf1, iabf2).real(), Wmat(iabf1, iabf2).imag()); // debug
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return Mmat;
}


matrix<cplxdb> compute_representation_on_equiv_k(const vec<double> &kprime, const matrix<cplxdb> &cmat_at_kprime,
                                                 const matrix<double> &lattice,
                                                 const matrix<double> &positions,
                                                 const vector<int> &atom_types,
                                                 const matrix<int> &rotmat_spg, const vec<double> &transi_spg,
                                                 const std::map<int, std::vector<abf_id>> &map_type_abfs_in,
                                                 const CODE_CHOICE & choice)
{
    auto wmat = compute_M_matrix(lattice, positions, atom_types, kprime, rotmat_spg, transi_spg, map_type_abfs_in, choice);
    return wmat * cmat_at_kprime * transpose(wmat, true);
}
