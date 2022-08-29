#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <string>
#include "spglib_utils.h"
#include "cell.h"
#include "io.h"
#include "kmesh.h"
#include "trans.h"

using std::cout;
using std::endl;
using std::logic_error;

int main (int argc, char *argv[])
{
    /* =============== */
    /* handling inputs */
    /* =============== */
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
        printf("euler angle: %8.5f %8.5f %8.5f %d\n", euler[0], euler[1], euler[2], is_proper);
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

    if (argc == 9)
    {
        vec<double> KRcoord(3);

        // one-line case
        code_choice = parse_code_choice(argv[3]);
        krmode = parse_krmode(argv[4]);
        printf("RSH choice: %s\n", CODE_CHOICE_STR[code_choice]);
        printf("      Mode: %s\n", KRMODE_STR[krmode]);

        for (int i = 0; i < 3; i++)
            KRcoord[i] = decode_fraction(argv[5+i]);
        mtxfile = argv[8];
        cout << "Reading matrix data from file: " << mtxfile << endl;
        const auto mat = read_mtx_cplxdb(mtxfile);
        printf("matrix size = (%d, %d), basis size = %d\n", mat.nr, mat.nc, basis.get_number_of_total_abfs());
        if (mat.nr != basis.get_number_of_total_abfs())
        {
            throw std::invalid_argument("The matrix size is not equal to the basis size, check input");
        }

        if (krmode == KRMODE::K)
        {
            printf("Computing all equivalent k-points to (%f, %f, %f)\n", KRcoord[0], KRcoord[1], KRcoord[2]);
            get_all_equiv_k(KRcoord, spgds.rotations, ks, isymops, code_choice);
            printf("Found %zu equivalent kpoints (including itself)\n", ks.size());
            // for (int ik = 0; ik < ks.size(); ik++)
            //     printf("%2d: %f %f %f\n", ik, ks[ik][0], ks[ik][1], ks[ik][2]);
            for (int ik = 0; ik < ks.size(); ik++)
            {
                printf("Transforming matrix at %7.4f %7.4f %7.4f -> %7.4f %7.4f %7.4f (%3d)\n",
                       KRcoord[0], KRcoord[1], KRcoord[2], ks[ik][0], ks[ik][1], ks[ik][2], ik);
                const auto euler = get_Euler_from_sym_matrix_spg(spgds.rotations[isymops[ik]], spgds.lattice, is_proper);
                if (is_proper)
                    printf("  Euler angle of Op %2d: %10.6f %10.6f %10.6f\n", isymops[ik]+1, euler[0], euler[1], euler[2]);
                else
                    printf("  Euler angle of Op %2d: %10.6f %10.6f %10.6f (plus inversion)\n", isymops[ik]+1, euler[0], euler[1], euler[2]);
                const auto k_mat = compute_representation_on_equiv_k(KRcoord, mat, spgds.lattice, spgds.positions, spgds.types,
                                                                     spgds.rotations[isymops[ik]], spgds.translations[isymops[ik]], map_type_abfs,
                                                                     code_choice);
                transformed_mat_outmtxfn = "abf_trans_out_equivk_" + std::to_string(ik) + "_symop_" + std::to_string(isymops[ik]+1) + ".mtx";
                sprintf(comment, "equiv k: %f %f %f", ks[ik][0], ks[ik][1], ks[ik][2]);
                write_mtx_cplxdb(k_mat, transformed_mat_outmtxfn, comment);
            }
            ks.clear();
            isymops.clear();
        }
    }

    if (argc == 4)
    {
        vector<vec<double>> krpoints;
        std::array<int, 3> ngs;
        vector<string> mtxfns;
        vector<matrix<cplxdb>> matrices;
        double output_thres = 1e-4; // negative to output anyway
        cout << "Reading code choice, K/R mode and matrice information from file: " << argv[3] << endl;
        read_matrix_inputs(argv[3], code_choice, krmode, ngs, krpoints, mtxfns, matrices);
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
                get_all_equiv_k(vk, spgds.rotations, ks, isymops, code_choice); // all ks that ~k = Vk
                for (int ik = 0; ik < ks.size(); ik++)
                {
                    const auto &k = ks[ik]; // k
                    if (!kgrids.have_irk(vk) or !kgrids.have_irk(k)) continue; // debug, only mapping IBZ kpts to itself
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
                    const auto &mat_k = matrices[ik_in_krpoints];
                    const auto &mat_vk = matrices[ivk];

                    const auto isymops = get_symops_connecting_k1_k2(k, vk, spgds.rotations);
                    cout << "All operations connecting two k: ";
                    for (const auto &iop: isymops)
                        cout << iop+1 << " ";
                    cout << endl;


                    for (auto isymop: isymops)
                    {
                        int ik_trans_from, ik_trans_to;
                        // consistency check
                        assert(is_same_k(vk, inverse(transpose(to_double(spgds.rotations[isymop]))) * k));
                        // do the transform, general
                        // ik_trans_from = ivk_in_grids;
                        // ik_trans_to = ik_in_grids;
                        // auto mmat = compute_M_matrix(spgds.lattice, spgds.positions, spgds.types, k,
                        //                              spgds.rotations[isymop], spgds.translations[isymop], map_type_abfs, code_choice);
                        // const auto mat_transformed = mmat * mat_vk * transpose(mmat, true);
                        // double maxabs_diff = maxabs(mat_transformed - mat_k);

                        // check the aims implementation particularly
                        // ik_trans_from = ik_in_grids;
                        // ik_trans_to = ivk_in_grids;
                        // auto mmat = compute_M_matrix_aims(spgds.lattice, spgds.positions, spgds.types, k,
                        //                                   spgds.rotations[isymop], spgds.translations[isymop], map_type_abfs);
                        // const auto mat_transformed = mmat * mat_k * transpose(mmat, true);
                        // double maxabs_diff = maxabs(mat_transformed - mat_vk);

                        // check the abacus implementation particularly
                        ik_trans_from = ik_in_grids;
                        ik_trans_to = ivk_in_grids;
                        auto mmat = compute_M_matrix_abacus(spgds.lattice, spgds.positions, spgds.types, k,
                                                            spgds.rotations[isymop], spgds.translations[isymop], map_type_abfs);
                        const auto mat_transformed = mmat * mat_k * transpose(mmat, true);
                        double maxabs_diff = maxabs(mat_transformed - mat_vk);

                        // print out the result
                        printf("    Sym. Op. %2d, |M(trans) - M(origi)|_max = %8.5f\n", isymop+1, maxabs_diff);
                        transformed_mat_outmtxfn = "abf_trans_out_ik_" + std::to_string(ik_trans_to+1) + 
                            "_from_" + std::to_string(ik_trans_from+1) + "_symop_" + std::to_string(isymop+1) + ".mtx";
                        mmat_outmtxfn = "M_ik_" + std::to_string(ik_trans_to+1) + 
                            "_from_" + std::to_string(ik_trans_from+1) + "_symop_" + std::to_string(isymop+1) + ".mtx";
                        write_mtx_cplxdb(mat_transformed, transformed_mat_outmtxfn);
                        write_mtx_cplxdb(mmat, mmat_outmtxfn);
                    }

                    // debug: LiF case, map IBZ k2 to BZ k7
                    // if (ik_in_grids != 1 || ik_equiv_in_grids != 6) continue;

                    // if (ik == ik_equiv_in_krpoints) // debug: mapping to itself
                        // const auto euler = get_Euler_from_sym_matrix_spg(spgds.rotations[isymop], spgds.lattice, is_proper);
                        /* wmat = transpose(wmat); // debug, test operation 23 */
                        /* wmat = transpose(wmat, true); // debug, test operation 23 */
                        // printf("Found k-point matrix mapping %4d(%6.3f, %6.3f %6.3f) -> %4d(%6.3f, %6.3f %6.3f) by symop. %2d, |Mk - Mkeq|_max = %f, |W Mk W^H - Mkeq|_max  = %f\n",
                        //        ik+1, kprime[0], kprime[1], kprime[2],
                        //        ik_equiv_in_krpoints+1, k_equiv[0], k_equiv[1], k_equiv[2],
                        //        isymop+1, maxabs(mat_kprime - mat_k_equiv), maxabs_diff);
                        // if (maxabs_diff > output_thres)
                            // debug
                            // print the symmetry operation
                            // cout << "     Symop " << isymop+1 << ": " << endl << spgdataset.get_operation_str_matform(isymop) << endl;
                            // const int i_invop = spgdataset.inverse_operation[isymop];
                            // cout << "Inverse op " << i_invop+1 << ": " << endl << spgdataset.get_operation_str_matform(i_invop) << endl;
                                /* wmat = transpose(wmat); // debug, test operation 23 */
                                /* wmat = transpose(wmat, true); // debug, test operation 23 */
                }
                isymops.clear();
                ks.clear();
            }
        }
    }


    return 0;
}
