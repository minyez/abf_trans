#include <algorithm>
#include <vector>
#include "kmesh.h"
#include "cell.h"

matrix<double> get_kgrids(std::array<int, 3> nks)
{
    const int nkpts = nks[0] * nks[1] * nks[2];
    std::vector<vec<double>> kgrids(nkpts);
    int ikpt = 0;
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

    return matrix<double>(kgrids);
}

void KGrids::set_kgrids()
{
    nkpts = nks[0] * nks[1] * nks[2];
    kpts = get_kgrids(nks);
}

KGrids::KGrids(int nkx, int nky, int nkz): nks({nkx, nky, nkz})
{
    set_kgrids();
}

KGrids::KGrids(const std::array<int, 3> &nks_in): nks(nks_in)
{
    set_kgrids();
}

void KGrids::rebuild_grids(int nkx, int nky, int nkz)
{
    nks = std::array<int, 3>{nkx, nky, nkz};
    set_kgrids();
    if (irkgrids_generated())
        clear_irkgrids();
}

void KGrids::rebuild_grids(const std::array<int, 3> &nks_in)
{
    nks = nks_in;
    set_kgrids();
    if (irkgrids_generated())
        clear_irkgrids();
}

void KGrids::clear_irkgrids()
{
    irkpts.resize(0, 0);
    indexmap_k_symop.clear();
    indexmap_k_irk.clear();
    irk_index.clear();
    irk_weights.clear();
    nirkpts = 0;
}

void KGrids::generate_irk_map(const SpgDS_c &dataset)
{
    vector<vec<double>> irkpts_vec;
    auto AAT = dataset.lattice * transpose(dataset.lattice);
    auto invAAT = inverse(AAT);

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
                auto tranmat = AAT * to_double(dataset.rotations[is]) * invAAT;
                auto vk = tranmat * kpt;
                move_to_center(vk, -0.5, false);
                // TODO: there might be precision problem
                const auto ite_vk = std::find(irkpts_vec.cbegin(), irkpts_vec.cend(), vk);
                if (ite_vk != irkpts_vec.cend())
                {
                    i_equiv_k = std::distance(irkpts_vec.cbegin(), ite_vk);
                    found_no_equiv = false;
                    break;
                }
            }
        if (found_no_equiv)
        {
            // a new irreducible k-points
            indexmap_k_irk.push_back(irkpts_vec.size());
            irkpts_vec.push_back(kpt);
            irk_index.push_back(ik);
            indexmap_k_symop.push_back(0);
            irk_weights.push_back(1);
        }
        else
        {
            indexmap_k_irk.push_back(i_equiv_k);
            indexmap_k_symop.push_back(is);
            irk_weights[i_equiv_k] += 1;
        }
    }

    nirkpts = irkpts_vec.size();
    irkpts = irkpts_vec;
}

void get_all_equiv_k(const vec<double> &k, const matrix<double> lattice,
                     const vector<matrix<int>> &rotmats_spg, vector<vec<double>> &equiv_ks, vector<int> &irots)
{
    if (!equiv_ks.empty()||!irots.empty())
    {
        std::cout << "Warning: clearing equiv_ks and irots" << std::endl;
        equiv_ks.clear();
        irots.clear();
    }
    auto AAT = lattice * transpose(lattice);
    auto invAAT = inverse(AAT);
    for (int isymop = 0; isymop < rotmats_spg.size(); isymop++)
    {
        const matrix<int> rotmat_spg = rotmats_spg[isymop];
        auto equiv_k = (AAT * to_double(rotmat_spg) * invAAT) * k;
        auto ite_ek = std::find(equiv_ks.cbegin(), equiv_ks.cend(), equiv_k);
        if (ite_ek != equiv_ks.cend())
        {
            equiv_ks.push_back(equiv_k);
            irots.push_back(isymop);
        }
    }
}
