#include <iostream>
#include <cstdlib>
#include <ctime>
#include <string>
#include "spglib_utils.h"
#include "cell.h"
#include "io.h"

using namespace std;

int main (int argc, char *argv[])
{
    /* =============== */
    /* handling inputs */
    /* =============== */
    std::cout << "Command: " << std::endl << "  ";
    for (int i = 0; i < argc; i++)
        std::cout << argv[i] << " ";
    std::cout << std::endl;
    if (argc < 3)
        throw logic_error("required cell and basis inputs are not parsed");

    /* cell input*/
    std::cout << "Reading cell information from " << argv[1] << std::endl;
    read_cell(argv[1], latt, posi_frac, types);
    std::cout << "  Lattice parameters: " << std::endl;
    for (int i = 0; i < 3; i++)
        printf("   a%1d:  %12.7f %12.7f %12.7f\n", i+1, latt(i, 0), latt(i, 1), latt(i, 2));
    std::cout << "  Fractional position and types of atoms. " << posi_frac.nr << " atoms in total: " << std::endl;
    for (int i = 0; i < posi_frac.nr; i++)
        printf("  %3d:  %12.7f %12.7f %12.7f %6d\n", i+1, posi_frac(i, 0), posi_frac(i, 1), posi_frac(i, 2), types[i]);
    /* post-processing atom type information */
    generate_map_type_iatom(types, inequiv_types, map_type_iatoms);
    std::cout << "  " << inequiv_types.size() << " inequivlent types: ";
    for (const auto &it: inequiv_types)
        std::cout << it << " ";
    std::cout << std::endl;

    /* Compute and display symmetry information */
    std::cout << "Compute symmetry related information by Spglib" << std::endl;
    SpgDS_c spgdataset(latt, posi_frac, types, 1.0e-5);
    spgdataset.show();

    /* basis input*/
    std::cout << std::endl << "Reading atomic basis information from " << argv[2] << std::endl;
    read_abf_ids(argv[2], inequiv_types, map_type_abfs);


    return 0;
}
