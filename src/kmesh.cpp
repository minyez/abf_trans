#include <algorithm>
#include <vector>
#include <iostream>
#include "kmesh.h"
#include "cell.h"
#include "mathtools.h"

matrix<double> get_kgrids(std::array<int, 3> nks, const CODE_CHOICE &code)
{
    const int nkpts = nks[0] * nks[1] * nks[2];
    std::vector<vec<double>> kgrids(nkpts);
    int ikpt = 0;
    switch (code)
    {
        case CODE_CHOICE::ORIG:
            for (int ikz = 0; ikz < nks[2]; ikz++)
            {
                double kz = double(ikz) / nks[2];
                if (kz > 0.5) kz -= 1.0;
                for (int iky = 0; iky < nks[1]; iky++)
                {
                    double ky = double(iky) / nks[1];
                    if (ky > 0.5) ky -= 1.0;
                    for (int ikx = 0; ikx < nks[0]; ikx++)
                    {
                        double kx = double(ikx) / nks[0];
                        if (kx > 0.5) kx -= 1.0;
                        kgrids[ikpt].resize(3);
                        kgrids[ikpt][0] = kx;
                        kgrids[ikpt][1] = ky;
                        kgrids[ikpt][2] = kz;
                        ikpt++;
                    }
                }
            }
            break;
        // by default aims does not have negative k-point.
        case CODE_CHOICE::AIMS:
            for (int ikx = 0; ikx < nks[0]; ikx++)
            {
                double kx = double(ikx) / nks[0];
                for (int iky = 0; iky < nks[1]; iky++)
                {
                    double ky = double(iky) / nks[1];
                    for (int ikz = 0; ikz < nks[2]; ikz++)
                    {
                        double kz = double(ikz) / nks[2];
                        kgrids[ikpt].resize(3);
                        kgrids[ikpt][0] = kx;
                        kgrids[ikpt][1] = ky;
                        kgrids[ikpt][2] = kz;
                        ikpt++;
                    }
                }
            }
            break;
    }

    return matrix<double>(kgrids);
}

void KGrids::set_kgrids()
{
    nkpts = nks[0] * nks[1] * nks[2];
    kpts = get_kgrids(nks, code);
}

KGrids::KGrids(int nkx, int nky, int nkz, const CODE_CHOICE &code_in): nks({nkx, nky, nkz}), code(code_in)
{
    set_kgrids();
}

KGrids::KGrids(const std::array<int, 3> &nks_in, const CODE_CHOICE &code_in): nks(nks_in), code(code_in)
{
    set_kgrids();
}

void KGrids::rebuild_grids(int nkx, int nky, int nkz, const CODE_CHOICE &code_in)
{
    nks = std::array<int, 3>{nkx, nky, nkz};
    code = code_in;
    set_kgrids();
    if (irkgrids_generated())
        clear_irkgrids();
}

void KGrids::rebuild_grids(const std::array<int, 3> &nks_in, const CODE_CHOICE &code_in)
{
    nks = nks_in;
    code = code_in;
    set_kgrids();
    if (irkgrids_generated())
        clear_irkgrids();
}

void KGrids::clear_irkgrids()
{
    irkpts.resize(0, 0);
    indexmap_k_symop.clear();
    indexmap_k_invsymop.clear();
    indexmap_k_irk.clear();
    irk_index.clear();
    irk_weights.clear();
    map_iirk_iks.clear();
    nirkpts = 0;
}

void KGrids::generate_irk_map(const SpgDS_c &dataset)
{
    vector<vec<double>> irkpts_vec;
    // auto AAT = dataset.lattice * transpose(dataset.lattice);
    // auto AAT = transpose(dataset.lattice) * dataset.lattice;
    // auto invAAT = inverse(AAT);

    if (irkgrids_generated())
        clear_irkgrids();

    for (int ik = 0; ik < kpts.nr; ik++)
    {
        auto kpt = kpts.get_row(ik);
        bool found_no_equiv = true;
        int i_equiv_k = -1;
        matrix<double> drot(3, 3);
        int is = 0;
        if (irkpts_vec.size() > 0)
            for (is = 0; is < dataset.n_operations; is++)
            {
                auto vk = inverse(transpose(to_double(dataset.rotations[is]))) * kpt;
                switch (code)
                {
                    case CODE_CHOICE::ORIG:
                        move_to_center(vk, -0.5, false); break;
                    case CODE_CHOICE::AIMS:
                        move_to_center(vk, 0.0, true); break;
                        /* move_to_center(vk, -0.5, false); break; */
                }
                for (int irk = 0; irk < irkpts_vec.size(); irk++)
                {
                    if (is_same_k(irkpts_vec[irk], vk))
                    {
                        i_equiv_k = irk;
                        found_no_equiv = false;
                        break;
                    }
                }
            }
        if (found_no_equiv)
        {
            // a new irreducible k-points
            indexmap_k_irk.push_back(irkpts_vec.size());
            irkpts_vec.push_back(kpt);
            irk_index.push_back(ik);
            indexmap_k_symop.push_back(0); // always the identity operation
            indexmap_k_invsymop.push_back(0); // always the identity operation
            irk_weights.push_back(1);
            map_iirk_iks[ik].push_back(ik);
        }
        else
        {
            indexmap_k_irk.push_back(i_equiv_k);
            indexmap_k_symop.push_back(is);
            indexmap_k_invsymop.push_back(dataset.inverse_operation[is]);
            irk_weights[i_equiv_k] += 1;
            map_iirk_iks[irk_index[i_equiv_k]].push_back(ik);
        }
    }

    // internal check
    int n_mapped_ik = 0;
    for (const auto &iirk: irk_index)
        n_mapped_ik += map_iirk_iks[iirk].size();
    assert(n_mapped_ik == nkpts);

    nirkpts = irkpts_vec.size();
    irkpts = irkpts_vec;
}


void get_all_equiv_k(const vec<double> &k, const matrix<double> lattice,
                     const vector<matrix<int>> &rotmats_spg, vector<vec<double>> &equiv_ks, vector<int> &irots,
                     const CODE_CHOICE &code)
{
    if (!equiv_ks.empty()||!irots.empty())
    {
        std::cout << "Warning: clearing equiv_ks and irots" << std::endl;
        equiv_ks.clear();
        irots.clear();
    }
    // const auto AAT = lattice * transpose(lattice);
    // const auto invAAT = inverse(AAT);
    for (int isymop = 0; isymop < rotmats_spg.size(); isymop++)
    {
        const matrix<int> &rotmat_spg = rotmats_spg[isymop];
        auto equiv_k = inverse(transpose(to_double(rotmat_spg))) * k;
        move_k_back(equiv_k, code);
        bool found = false;
        for (const auto &already_found_ek: equiv_ks)
        {
            if (is_same_k(already_found_ek, equiv_k))
            {
                found = true;
                break;
            }
        }
        if (!found)
        {
            equiv_ks.push_back(equiv_k);
            irots.push_back(isymop);
        }
    }
}

vector<int> get_all_symops_connecting_ks(const vec<double> &k, const vec<double> &Vk,
                                         const matrix<double> lattice,
                                         const vector<matrix<int>> &rotmats_spg)
{
    vector<int> indices;
    // const auto AAT = lattice * transpose(lattice);
    // const auto invAAT = inverse(AAT);
    for (int isymop = 0; isymop < rotmats_spg.size(); isymop++)
    {
        const matrix<int> &rotmat_spg = rotmats_spg[isymop];
        auto equiv_k = inverse(transpose(to_double(rotmat_spg))) * k;
        // std::cout << isymop << " " << k << " " << equiv_k << " " << Vk << std::endl; // debug
        if (is_same_k(equiv_k, Vk))
        {
            indices.push_back(isymop);
        }
    }
    return indices;
}

matrix<double> move_k_back(matrix<double> &kpts, const CODE_CHOICE &code)
{
    matrix<double> K(kpts.nr, kpts.nc);
    switch (code)
    {
        case CODE_CHOICE::ORIG:
            for (int ia = 0; ia < K.nr; ia++)
                for (int ic = 0; ic < K.nc; ic++)
                {
                    K(ia, ic) = shift_to_unit(kpts(ia, ic), -0.5, false);
                }
            break;
        case CODE_CHOICE::AIMS:
            for (int ia = 0; ia < K.nr; ia++)
                for (int ic = 0; ic < K.nc; ic++)
                {
                    K(ia, ic) = shift_to_unit(kpts(ia, ic), 0.0, true);
                }
            break;
    }
    return K;
}

vec<double> move_k_back(vec<double> &kpts, const CODE_CHOICE &code)
{
    vec<double> K(kpts.n);
    switch (code)
    {
        case CODE_CHOICE::ORIG:
                for (int ic = 0; ic < K.n; ic++)
                    K[ic] = shift_to_unit(kpts[ic], -0.5, false);
            break;
        case CODE_CHOICE::AIMS:
                for (int ic = 0; ic < K.n; ic++)
                    K[ic] = shift_to_unit(kpts[ic], 0.0, true);
            break;
    }
    return K;
}
