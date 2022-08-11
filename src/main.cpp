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
    SpgDS_c spgdataset(latt, posi_frac, types, 1.0e-5);
    spgdataset.show(true);
    /* cout << "Standard positions in original lattice:" << endl; */
    /* cout << (spgdataset.std_positions - spgdataset.origin_shift )* inverse(spgdataset.transformation_matrix); */
    /* cout << spgdataset.std_positions * inverse(spgdataset.transformation_matrix); */

    /* basis input*/
    cout << endl << "Reading atomic basis information from " << argv[2] << endl;
    read_abf_ids(argv[2], inequiv_types, map_type_abfs);
    ABF basis(types, map_type_abfs);

    /* prepare required quantities*/
    vector<vec<double>> equiv_ks;
    vector<int> isymops;
    bool is_proper;
    char comment[80];
    string mtxfile, outmtxfn;

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
            get_all_equiv_k(KRcoord, spgdataset.lattice, spgdataset.rotations, equiv_ks, isymops, code_choice);
            printf("Found %zu equivalent kpoints (including itself)\n", equiv_ks.size());
            // for (int ik = 0; ik < ks.size(); ik++)
            //     printf("%2d: %f %f %f\n", ik, ks[ik][0], ks[ik][1], ks[ik][2]);
            for (int ik = 0; ik < equiv_ks.size(); ik++)
            {
                printf("Transforming matrix at %7.4f %7.4f %7.4f -> %7.4f %7.4f %7.4f (%3d)\n",
                       KRcoord[0], KRcoord[1], KRcoord[2], equiv_ks[ik][0], equiv_ks[ik][1], equiv_ks[ik][2], ik);
                const auto euler = get_Euler_from_sym_matrix_spg(spgdataset.rotations[isymops[ik]], spgdataset.lattice, is_proper);
                if (is_proper)
                    printf("  Euler angle of Op %2d: %10.6f %10.6f %10.6f\n", isymops[ik]+1, euler[0], euler[1], euler[2]);
                else
                    printf("  Euler angle of Op %2d: %10.6f %10.6f %10.6f (plus inversion)\n", isymops[ik]+1, euler[0], euler[1], euler[2]);
                const auto k_mat = compute_representation_on_equiv_k(KRcoord, mat, spgdataset.lattice, spgdataset.positions, spgdataset.types,
                                                                     spgdataset.rotations[isymops[ik]], spgdataset.translations[isymops[ik]], map_type_abfs,
                                                                     code_choice);
                outmtxfn = "abf_trans_out_equivk_" + std::to_string(ik) + "_symop_" + std::to_string(isymops[ik]+1) + ".mtx";
                sprintf(comment, "equiv k: %f %f %f", equiv_ks[ik][0], equiv_ks[ik][1], equiv_ks[ik][2]);
                write_mtx_cplxdb(k_mat, outmtxfn, comment);
            }
            equiv_ks.clear();
            isymops.clear();
        }
    }
    else if (argc == 4)
    {
        vector<vec<double>> krpoints;
        std::array<int, 3> ngs;
        vector<string> mtxfns;
        vector<matrix<cplxdb>> matrices;
        double output_thres = 1e-4; // negative to output anyway
        cout << "Reading code choice, K/R mode and matrice information from file: " << argv[3] << endl;
        read_matrix_inputs(argv[3], ngs, krpoints, mtxfns, matrices);
        if (krmode == KRMODE::K)
        {
            KGrids kgrids(ngs, code_choice);
            kgrids.generate_irk_map(spgdataset);
            for (int ik = 0; ik < krpoints.size(); ik++)
            {
                const auto &kprime = krpoints[ik]; // k'
                const auto ik_in_grids = kgrids.index(kprime);
                // only the IBZ kpts
                // if (!kgrids.have_irk(kprime)) continue;
                get_all_equiv_k(kprime, spgdataset.lattice, spgdataset.rotations, equiv_ks, isymops, code_choice); // all ks that k=Vk'
                for (int ik_equiv = 0; ik_equiv < equiv_ks.size(); ik_equiv++)
                {
                    const auto &k_equiv = equiv_ks[ik_equiv]; // k
                    const auto &isymop = isymops[ik_equiv];

                    const auto ite_k_equiv_in_krpoints = std::find(krpoints.cbegin(), krpoints.cend(), equiv_ks[ik_equiv]);
                    if (ite_k_equiv_in_krpoints == krpoints.cend()) continue;
                    const auto &ik_equiv_in_grids = kgrids.index(k_equiv);
                    int ik_equiv_in_krpoints = std::distance(krpoints.cbegin(), ite_k_equiv_in_krpoints);
                    if (ik != ik_equiv_in_krpoints)
                    {
                        const auto euler = get_Euler_from_sym_matrix_spg(spgdataset.rotations[isymop], spgdataset.lattice, is_proper);
                        const auto &mat_kprime = matrices[ik];
                        const auto &mat_k_equiv = matrices[ik_equiv_in_krpoints];
                        const auto wmat = compute_W_matrix(spgdataset.lattice, spgdataset.positions, spgdataset.types, kprime,
                                                           spgdataset.rotations[isymop],
                                                           spgdataset.translations[isymop], map_type_abfs, code_choice);
                        const auto mat_k_transformed = wmat * mat_kprime * transpose(wmat, true);
                        double maxabs_diff = maxabs(mat_k_transformed - mat_k_equiv);
                        // printf("Found k-point matrix mapping %4d(%6.3f, %6.3f %6.3f) -> %4d(%6.3f, %6.3f %6.3f) by symop. %2d, |Mk - Mkeq|_max = %f, |W Mk W^H - Mkeq|_max  = %f\n",
                        //        ik+1, kprime[0], kprime[1], kprime[2],
                        //        ik_equiv_in_krpoints+1, k_equiv[0], k_equiv[1], k_equiv[2],
                        //        isymop+1, maxabs(mat_kprime - mat_k_equiv), maxabs_diff);
                        printf("Found k-point matrix mapping %4d(%6.3f, %6.3f %6.3f) -> %4d(%6.3f, %6.3f %6.3f)\n",
                               ik_in_grids+1, kprime[0], kprime[1], kprime[2],
                               ik_equiv_in_grids+1, k_equiv[0], k_equiv[1], k_equiv[2]);
                        printf("    by symop. %2d, |Mk - Mkeq|_max = %8.5e, |W Mk W^H - Mkeq|_max  = %8.5e",
                               isymop+1, maxabs(mat_kprime - mat_k_equiv), maxabs_diff);
                        const auto iops_kpke = get_all_symops_connecting_ks(kprime, k_equiv, spgdataset.lattice, spgdataset.rotations);
                        if (maxabs_diff > output_thres)
                            printf(" (X)");
                        printf("\n");
                        cout << "All operations connecting two k: ";
                        for (const auto &iop: iops_kpke)
                            cout << iop+1 << " ";
                        cout << endl;
                        if (maxabs_diff > output_thres)
                        {
                            // debug
                            outmtxfn = "abf_trans_out_ik_" + std::to_string(ik_equiv_in_grids+1) + 
                                "_from_" + std::to_string(ik_in_grids+1) + "_symop_" + std::to_string(isymop+1) + ".mtx";
                            write_mtx_cplxdb(mat_k_transformed, outmtxfn);
                            // print the symmetry operation
                            // cout << "     Symop " << isymop+1 << ": " << endl << spgdataset.get_operation_str_matform(isymop) << endl;
                            // const int i_invop = spgdataset.inverse_operation[isymop];
                            // cout << "Inverse op " << i_invop+1 << ": " << endl << spgdataset.get_operation_str_matform(i_invop) << endl;
                            for (const auto &iop: iops_kpke)
                            {
                                if (iop == isymop) continue;
                                const auto wmat = compute_W_matrix(spgdataset.lattice, spgdataset.positions, spgdataset.types, kprime,
                                                                   spgdataset.rotations[iop],
                                                                   spgdataset.translations[iop], map_type_abfs, code_choice);
                                const auto mat_k_transformed = wmat * mat_kprime * transpose(wmat, true);
                                double maxabs_diff = maxabs(mat_k_transformed - mat_k_equiv);
                                printf("    by symop. %2d,                                |W Mk W^H - Mkeq|_max  = %8.5e", iop+1, maxabs_diff);
                                if (maxabs_diff > output_thres)
                                    printf(" (X)");
                                printf("\n");
                            }
                        }
                    }
                }
                isymops.clear();
                equiv_ks.clear();
            }
        }
    }


    return 0;
}
