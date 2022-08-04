#include <iostream>
#include <cstdlib>
#include <ctime>
#include <string>
#include "spglib_utils.h"
#include "cell.h"
#include "io.h"

using std::cout;
using std::endl;
using std::logic_error;

int main (int argc, char *argv[])
{
    /* =============== */
    /* handling inputs */
    /* =============== */
    cout << "Command: " << endl << "  ";
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

    if (argc == 8)
    {
        vec<double> kprime(3);
        for (int i = 0; i < 3; i++)
            kprime[i] = std::stod(argv[5+i]);
    }


    return 0;
}
