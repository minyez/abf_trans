#include <iostream>
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
    cout << "Command: " << argc << endl << "  ";
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
    spgdataset.show();
    /* cout << "Standard positions in original lattice:" << endl; */
    /* cout << (spgdataset.std_positions - spgdataset.origin_shift )* inverse(spgdataset.transformation_matrix); */
    /* cout << spgdataset.std_positions * inverse(spgdataset.transformation_matrix); */

    /* basis input*/
    cout << endl << "Reading atomic basis information from " << argv[2] << endl;
    read_abf_ids(argv[2], inequiv_types, map_type_abfs);
    ABF basis(types, map_type_abfs);

    /* prepare required quantities*/
    string choice, mode, mtxfile;
    vec<double> KRcoord(3);
    vector<vec<double>> ks;
    vector<int> isymops;
    string outmtxfn;
    char comment[80];

    if (argc == 9)
    {
        // one-line case
        choice = argv[3];
        mode = argv[4];
        printf("RSH choice: %s\n", choice.c_str());
        printf("      Mode: %s\n", mode.c_str());

        for (int i = 0; i < 3; i++)
            KRcoord[i] = decode_fraction(argv[5+i]);
        mtxfile = argv[8];
        const auto mat = read_mtx_cplxdb(mtxfile);

        if (mode == "K")
        {
            printf("Computing all equivalent k-points to (%f, %f, %f)\n", KRcoord[0], KRcoord[1], KRcoord[2]);
            get_all_equiv_k(KRcoord, spgdataset.lattice, spgdataset.rotations, ks, isymops);
            printf("Found %zu equivalent kpoints (including itself)\n", ks.size());
            for (int ik = 0; ik < ks.size(); ik++)
            {
                const auto k_mat = compute_representation_on_equiv_k(KRcoord, mat, spgdataset.lattice, spgdataset.positions, spgdataset.types, 
                                                                     spgdataset.rotations[isymops[ik]], spgdataset.translations[isymops[ik]], map_type_abfs);
                outmtxfn = "abf_trans_out_equivk_" + std::to_string(ik) + "_symop_" + std::to_string(isymops[ik]) + ".mtx";
                sprintf(comment, "equiv k: %f %f %f", ks[ik][0], ks[ik][1], ks[ik][2]);
                write_mtx_cplxdb(k_mat, outmtxfn);
            }
            ks.clear();
            isymops.clear();
        }
    }


    return 0;
}
