#include "../src/spglib_utils.h"
#include <iostream>
#include <string>

using namespace std;

void SpgDS_c_internal_consistency_check(const string & msg, const SpgDS_c &dataset)
{
    cout << "Checking " << msg << "..." << endl;
    auto latt_from_std = dataset.transformation_matrix * dataset.std_lattice;
    cout << "Current lattice from input:" << endl;
    cout << dataset.lattice;
    cout << "Current lattice by transformation_matrix * std_lattice:" << endl;
    cout << latt_from_std;
    assert(dataset.lattice == latt_from_std);
    cout << "Idealized lattice:" << endl;
    cout << dataset.get_ideal_lattice();
    cout << "Number of atoms in idealized standard lattice = " << dataset.n_std_atoms << endl;
    cout << "Atom positions :" << endl;
    cout << dataset.std_positions;
    cout << "Convert atom positions in input cell to ideal lattice:" << endl;
    cout << dataset.positions * dataset.transformation_matrix + dataset.origin_shift;
}

void check_prim_diamond()
{
    const double halfa = 1.78356;
    double latt[9] = { 0, halfa, halfa, halfa, 0, halfa, halfa, halfa, 0};
    double posi[6] = {0, 0, 0, 0.25, 0.25, 0.25};
    SpgDS_c dataset({3, 3, latt}, {2, 3, posi}, {6, 6});
    SpgDS_c_internal_consistency_check("Primitive C diamond", dataset);
}

void check_nonideal_prim_diamond()
{
    const double halfa = 1.78356;
    const double a = 2 * halfa;
    double latt[9] = { 0, halfa, halfa, halfa, halfa, 2*halfa, halfa, halfa, 0};
    double posi[6] = {0, 0, 0, 0.0, 0.25, 0.25};
    SpgDS_c dataset({3, 3, latt}, {2, 3, posi}, {6, 6});
    SpgDS_c_internal_consistency_check("Non-ideal Primitive diamond", dataset);

    double ideal_latt[9] = {a, 0, 0, 0, a, 0, 0, 0, a};
    matrix<double> ideal_lattm(3, 3, ideal_latt);
}

int main (int argc, char *argv[])
{
    check_prim_diamond();
    check_nonideal_prim_diamond();
    return 0;
}
