#include "../src/kmesh.h"
#include "testutils.h"
#include <iostream>

using namespace std;

void test_get_kgrids()
{
    matrix<double> kgrids;

    cout << "kgrids for (2, 2, 2)" << endl;
    kgrids = get_kgrids({2, 2, 2});
    assert(kgrids.nr == 8);
    // cout << kgrids;

    cout << "kgrids for (3, 3, 3)" << endl;
    kgrids = get_kgrids({3, 3, 3});
    assert(kgrids.nr == 27);
    // cout << kgrids;

    cout << "kgrids for (4, 2, 3)" << endl;
    kgrids = get_kgrids({4, 2, 3});
    assert(kgrids.nr == 24);
    // cout << kgrids;

    cout << "kgrids for (5, 5, 5)" << endl;
    kgrids = get_kgrids({5, 5, 5});
    assert(kgrids.nr == 125);
    // cout << kgrids;
}

void test_KGrids_irkgrids_diamond_orig()
{
    // diamond 2x2x2
    matrix<double> latt(3, 3), posi(2, 3);
    latt(0, 1) = latt(0, 2) = latt(1, 0) = latt(1, 2) = latt(2, 0) = latt(2, 1) = 0.5;
    posi(1, 0) = posi(1, 1) = posi(1, 2) = 0.25;
    SpgDS_c dataset(latt, posi, {6, 6});
    KGrids kgrids(2, 2, 2, CODE_CHOICE::ORIG);
    kgrids.generate_irk_map(dataset);
    assert(kgrids.nirkpts == 3);
    cout << "Irreducible points for 2x2x2" << endl;
    for (int i = 0; i < kgrids.nirkpts; i++)
        cout << kgrids.irkpts.get_row(i) << " " << kgrids.irk_weights[i] << endl;

    kgrids.rebuild_grids(4, 2, 3, CODE_CHOICE::ORIG);
    kgrids.generate_irk_map(dataset);
    assert(kgrids.nirkpts == 12);
    cout << "Irreducible points for 4x2x3" << endl;
    for (int i = 0; i < kgrids.nirkpts; i++)
        cout << kgrids.irkpts.get_row(i) << " " << kgrids.irk_weights[i] << endl;
}

void test_KGrids_irkgrids_nacl_aims()
{
    const double onethird = 1./3.;
    const double twothird = 2./3.;
    const double kpts_k333_aims[] = {0.00000,  0.00000,  0.00000,
                                     0.00000,  0.00000,  onethird,
                                     0.00000,  0.00000,  twothird,
                                     0.00000,  onethird,  0.00000,
                                     0.00000,  onethird,  onethird,
                                     0.00000,  onethird,  twothird,
                                     0.00000,  twothird,  0.00000,
                                     0.00000,  twothird,  onethird,
                                     0.00000,  twothird,  twothird,
                                     onethird,  0.00000,  0.00000,
                                     onethird,  0.00000,  onethird,
                                     onethird,  0.00000,  twothird,
                                     onethird,  onethird,  0.00000,
                                     onethird,  onethird,  onethird,
                                     onethird,  onethird,  twothird,
                                     onethird,  twothird,  0.00000,
                                     onethird,  twothird,  onethird,
                                     onethird,  twothird,  twothird,
                                     twothird,  0.00000,  0.00000,
                                     twothird,  0.00000,  onethird,
                                     twothird,  0.00000,  twothird,
                                     twothird,  onethird,  0.00000,
                                     twothird,  onethird,  onethird,
                                     twothird,  onethird,  twothird,
                                     twothird,  twothird,  0.00000,
                                     twothird,  twothird,  onethird,
                                     twothird,  twothird,  twothird};
    matrix<double> kpts(27, 3, kpts_k333_aims);
    KGrids kgrids(3, 3, 3, CODE_CHOICE::AIMS);
    assert(kpts == kgrids.kpts);

    matrix<double> latt(3, 3), posi(2, 3);
    latt(0, 1) = latt(0, 2) = latt(1, 0) = latt(1, 2) = latt(2, 0) = latt(2, 1) = 0.5;
    posi(1, 0) = posi(1, 1) = posi(1, 2) = 0.50;
    vector<int> atom_types{11, 17};

    SpgDS_c dataset(latt, posi, atom_types);
    kgrids.generate_irk_map(dataset);
    const double irkpts_k333_aims[] = {0.00000,  0.00000,  0.00000,
                                       0.00000,  0.00000,  onethird,
                                       0.00000,  onethird,  onethird,
                                       0.00000,  onethird,  twothird};
    matrix<double> irkpts(4, 3, irkpts_k333_aims);
    assert(irkpts == kgrids.irkpts);
    cout << irkpts << kgrids.irkpts;

    printf("Checking irk indices ...\n");
    const vector<int> i_irks {0, 1, 4, 5};
    assert( vec<int>(i_irks) == vec<int>(kgrids.irk_index));
    
    printf("Checking iirk to ik mapping ...\n");
    map<int, vector<int>> map_iirk_iks;
    map_iirk_iks[0] = vector<int>({0});
    map_iirk_iks[1] = vector<int>({1, 2, 3, 6, 9, 13, 18, 26});
    map_iirk_iks[4] = vector<int>({4, 8, 10, 12, 20, 24});
    map_iirk_iks[5] = vector<int>({5, 7, 11, 14, 15, 16, 17, 19, 21, 22, 23, 25});
    for (auto iirk: kgrids.irk_index)
    {
        assert(vec<int>(map_iirk_iks[iirk]) == vec<int>(kgrids.map_iirk_iks[iirk]));
    }
}

void test_get_all_equiv_k_nacl()
{
    matrix<double> latt(3, 3), posi(2, 3);
    latt(0, 1) = latt(0, 2) = latt(1, 0) = latt(1, 2) = latt(2, 0) = latt(2, 1) = 0.5;
    posi(1, 0) = posi(1, 1) = posi(1, 2) = 0.50;
    vector<int> atom_types{11, 17};

    SpgDS_c dataset(latt, posi, atom_types);
    vec<double> kvec(3), ktarget(3);
    // 4( 0.000,  0.500,  0.500) ->    7( 0.500,  0.500,  0.000)
    /* kvec[0] = 0.0, kvec[1] = kvec[2] = 0.5; ktarget[2] = 0.0, ktarget[0] = ktarget[1] = 0.5; */
    // 7( 0.500,  0.500,  0.000) ->    4( 0.000,  0.500,  0.500)
    /* kvec[2] = 0.0, kvec[0] = kvec[1] = 0.5; kref[0] = 0.0, kref[1] = kref[2] = 0.5; */
    /* 6( 0.500,  0.000,  0.500) ->    4( 0.000,  0.500,  0.500) */
    // kvec[1] = 0.0, kvec[0] = kvec[2] = 0.5; ktarget[0] = 0.0, ktarget[1] = ktarget[2] = 0.5;
    /* 6( 0.500,  0.000,  0.500) ->    7( 0.500,  0.500,  0.000) */
    kvec[1] = 0.0, kvec[0] = kvec[2] = 0.5; ktarget[2] = 0.0, ktarget[1] = ktarget[0] = 0.5;
    vector<vec<double>> ks;
    vector<int> isymops;
    const matrix<double> AAT = dataset.lattice * transpose(dataset.lattice);
    const auto invAAT= inverse(AAT);
    cout << AAT << invAAT;
    get_all_equiv_k(kvec, latt, dataset.rotations, ks, isymops, CODE_CHOICE::AIMS);
    for (int ik_equiv = 0; ik_equiv < ks.size(); ik_equiv++)
    {
        const vec<double> &kequiv = ks[ik_equiv];
        if (!(kequiv == ktarget)) continue;
        const auto &isymop = isymops[ik_equiv];
        cout << "k: " << kvec << " -> " << "equivk: " << kequiv << " by symop " << isymop+1 << endl;
        // cout << isymop << endl;
        const auto a_v_ia = AAT * to_double(dataset.rotations[isymop]) * invAAT;
        cout << "V:\n" << dataset.rotations[isymop];
        cout << "A V A^-1:\n" << a_v_ia;
        const auto a_v_ia_k = a_v_ia * kvec;
        const auto a_iv_ia = AAT * to_double(dataset.rotations[dataset.inverse_operation[isymop]]) * invAAT;
        cout << "A V^-1 A^-1:\n" << a_iv_ia;
        const auto a_iv_ia_kequiv = a_iv_ia * kequiv;
        assert(vec_equal(kequiv, a_v_ia_k, 1e-5));
        assert(vec_equal(kvec, a_iv_ia_kequiv, 1e-5));
    }
}

int main (int argc, char *argv[])
{
    // test_get_kgrids();
    // test_KGrids_irkgrids_diamond();
    test_KGrids_irkgrids_nacl_aims();
    test_get_all_equiv_k_nacl();
    return 0;
}
