#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <string>
#include "cell.h"
#include "io.h"
#include "kmesh.h"
#include "logger.h"
#include "spglib_utils.h"
#include "trans.h"

using std::cout;
using std::endl;
using std::logic_error;

int main (int argc, char *argv[])
{
    const bool debug = false;
    /* =============== */
    /* handling inputs */
    /* =============== */
    logger.open("log.txt");

    cout << "Command: " << argc << " arguments" << endl << "  ";
    for (int i = 0; i < argc; i++)
        cout << argv[i] << " ";
    cout << endl;
    if (argc < 3)
        throw logic_error("required cell and basis inputs are not parsed");

    /* cell input*/
    cout << "Reading cell information from " << argv[1] << endl;
    read_cell(argv[1], latt, posi_frac, types);
    cout << "  Lattice parameters: " << endl;
    for (int i = 0; i < 3; i++)
        printf("   a%1d:  %12.7f %12.7f %12.7f\n", i+1, latt(i, 0), latt(i, 1), latt(i, 2));
    cout << "  Fractional position and types of atoms. " << posi_frac.nr << " atoms in total: " << endl;
    for (int i = 0; i < posi_frac.nr; i++)
        printf("  %3d:  %12.7f %12.7f %12.7f %6d\n", i+1, posi_frac(i, 0), posi_frac(i, 1), posi_frac(i, 2), types[i]);
    /* post-processing atom type information */
    generate_map_type_iatom(types, inequiv_types, map_type_iatoms);
    cout << "  " << inequiv_types.size() << " inequivlent types: ";
    for (const auto &it: inequiv_types)
        cout << it << " ";
    cout << endl;

    /* Compute and display symmetry information */
    cout << "Compute symmetry related information by Spglib" << endl;
    SpgDS_c spgds(latt, posi_frac, types, 1.0e-5);
    spgds.show(true, true);
    for (int iop = 0; iop < spgds.n_operations; iop++)
    {
        bool is_proper;
        auto euler = get_Euler_from_sym_matrix_spg(spgds.rotations[iop], spgds.lattice, is_proper);
        logger << "euler angle: " << euler[0] << " " << euler[1] << " " << euler[2] << " " << is_proper << endl;
    }
    /* cout << "Standard positions in original lattice:" << endl; */
    /* cout << (spgdataset.std_positions - spgdataset.origin_shift )* inverse(spgdataset.transformation_matrix); */
    /* cout << spgdataset.std_positions * inverse(spgdataset.transformation_matrix); */

    /* basis input*/
    cout << endl << "Reading atomic basis information from " << argv[2] << endl;
    read_abf_ids(argv[2], inequiv_types, map_type_abfs);
    ABF basis(types, map_type_abfs);
    cout << "Size of atomic basis: " << basis.get_number_of_total_abfs() << endl;
    /* prepare required quantities*/
    vector<vec<double>> ks;
    vector<int> isymops;
    bool is_proper;
    char comment[80];
    string mtxfile, transformed_mat_outmtxfn, mmat_outmtxfn;

    if (argc == 4)
    {
        vector<vec<double>> krpoints;
        std::array<int, 3> ngs;
        vector<string> matfns;
        vector<matrix<cplxdb>> matrices;
        cout << "Reading code choice, K/R mode and matrice information from file: " << argv[3] << endl;
        read_matrix_inputs(argv[3], code_choice, krmode, ngs, krpoints, matfns, matrices);
        if (krmode == KRMODE::K)
        {
            KGrids kgrids(ngs, code_choice);
            kgrids.generate_irk_map(spgds);
            for (int ivk = 0; ivk < krpoints.size(); ivk++)
            {
                const auto &vk = krpoints[ivk]; // Vk
                const auto ivk_in_grids = kgrids.index(vk);
                if (ivk_in_grids < 0) continue; // not found in the kgrids
                // if (ivk_in_grids != 0) continue; // debug, only the Gamma point
                // if (!kgrids.have_irk(vk)) continue; // debug, only IBZ kpts
                get_all_equiv_k(vk, spgds.rotations, ks, isymops, code_choice); // decompose Vk = V k
                for (int ik = 0; ik < ks.size(); ik++)
                {
                    const auto &k = ks[ik]; // k
                    // if (!kgrids.have_irk(vk) or !kgrids.have_irk(k)) continue; // debug, only mapping IBZ kpts to itself
                    const auto &ik_in_grids = kgrids.index(k);

                    int ik_in_krpoints;
                    for (ik_in_krpoints = 0; ik_in_krpoints < krpoints.size(); ik_in_krpoints++)
                    {
                        if (is_same_k(k, krpoints[ik_in_krpoints]))
                            break;
                    }
                    if (ik_in_krpoints == krpoints.size()) continue; // not found the transformed point in the grid

                    printf("Found k-point matrix mapping: k %4d(%6.3f, %6.3f %6.3f) -> Vk %4d(%6.3f, %6.3f %6.3f)\n",
                           ik_in_grids+1, k[0], k[1], k[2],
                           ivk_in_grids+1, vk[0], vk[1], vk[2]);
                    logger << "===== k " << ik_in_grids+1 << "(" << k << ") -> Vk " << ivk_in_grids+1 << "(" << vk << ") =====" << endl;
                    const auto &mat_k = matrices[ik_in_krpoints];
                    const auto &mat_vk = matrices[ivk];

                    const auto isymops = get_symops_connecting_k1_k2(k, vk, spgds.rotations);
                    logger << "All operations connecting two k: ";
                    for (const auto &iop: isymops)
                        logger << iop+1 << " ";
                    logger << endl;


                    for (auto isymop: isymops)
                    {
                        logger << "----- Sym. Op. " << isymop+1 << endl;
                        int ik_trans_from, ik_trans_to;
                        matrix<cplxdb> mmat, mat_transformed;
                        double maxabs_diff;
                        // consistency check
                        assert(is_same_k(vk, inverse(transpose(to_double(spgds.rotations[isymop]))) * k));

                        // do the transform
                        ik_trans_from = ik_in_grids;
                        ik_trans_to = ivk_in_grids;
                        mmat = compute_M_matrix(spgds.lattice, spgds.positions, spgds.types, k,
                                                spgds.rotations[isymop], spgds.translations[isymop], map_type_abfs, code_choice);
                        mat_transformed = mmat * mat_k * transpose(mmat, true);
                        maxabs_diff = maxabs(mat_transformed - mat_vk);

                        // print out the result
                        printf("    Sym. Op. %2d, |M(trans) - M(origi)|_max = %8.5f\n", isymop+1, maxabs_diff);
                        if (debug)
                        {
                            transformed_mat_outmtxfn = "abf_trans_out_ik_" + std::to_string(ik_trans_to+1) + 
                                "_from_" + std::to_string(ik_trans_from+1) + "_symop_" + std::to_string(isymop+1) + ".mtx";
                            mmat_outmtxfn = "M_ik_" + std::to_string(ik_trans_to+1) + 
                                "_from_" + std::to_string(ik_trans_from+1) + "_symop_" + std::to_string(isymop+1) + ".mtx";
                            write_mtx_cplxdb(mat_transformed, transformed_mat_outmtxfn);
                            write_mtx_cplxdb(mmat, mmat_outmtxfn);
                        }
                    }
                }
                isymops.clear();
                ks.clear();
            }
        }
    }

    logger.close();

    return 0;
}
