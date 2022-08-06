#include "../src/trans.h"
#include "../src/spglib_utils.h"
#include "testutils.h"
#include <iostream>

using namespace std;

void test_prim_nacl_spd()
{
    matrix<double> latt(3, 3), posi(2, 3);
    latt(0, 1) = latt(0, 2) = latt(1, 0) = latt(1, 2) = latt(2, 0) = latt(2, 1) = 0.5;
    posi(1, 0) = posi(1, 1) = posi(1, 2) = 0.50;
    vector<int> atom_types{11, 17};
    bool is_proper;

    SpgDS_c dataset(latt, posi, atom_types);
    dataset.show_cell();
    dataset.show();

    map_type_abfs[11].push_back({0});
    map_type_abfs[11].push_back({1});
    map_type_abfs[11].push_back({2});
    map_type_abfs[17].push_back({0});
    map_type_abfs[17].push_back({1});
    map_type_abfs[17].push_back({2});

    vec<double> kprime(3);
    const int i_symop = 2;
    kprime[0] = 0.5;
    cout << "k': "<< kprime << endl;
    cout << "rotation mat: "<< endl << dataset.rotations[i_symop];
    cout << "inverse: " << endl << inverse(to_double(dataset.rotations[i_symop]));
    cout << "transition v: " << dataset.translations[i_symop] << endl;
    auto euler = get_Euler_from_sym_matrix_spg(dataset.rotations[i_symop], latt, is_proper);
    printf("Euler angle: %f %f %f\n", euler[0], euler[1], euler[2]);
    auto Wmat = compute_W_matrix(latt, posi, atom_types, kprime,
                                 dataset.rotations[i_symop],
                                 dataset.translations[i_symop], map_type_abfs, "orig");
    // std::cout << Wmat;

    cout << "Testing unitary property " << endl;
    auto iden = Wmat * transpose(Wmat, true);
    matrix<cplxdb> iden_ref(iden.nr, iden.nc);
    iden_ref.set_diag(1);
    assert(iden == iden_ref);

    map_type_abfs.clear();
}

void test_prim_diamond_spd()
{
    matrix<double> latt(3, 3), posi(2, 3);
    latt(0, 1) = latt(0, 2) = latt(1, 0) = latt(1, 2) = latt(2, 0) = latt(2, 1) = 0.5;
    posi(1, 0) = posi(1, 1) = posi(1, 2) = 0.25;
    vector<int> atom_types{6, 6};
    bool is_proper;
    SpgDS_c dataset(latt, posi, atom_types);
    dataset.show_cell();
    dataset.show();

    map_type_abfs[6].push_back({0});
    map_type_abfs[6].push_back({1});
    map_type_abfs[6].push_back({2});

    vec<double> kprime(3);
    // FIXME: cannot map to a correct Mu'
    const int i_symop = 3;
    kprime[0] = 0.5;
    cout << "k': "<< kprime << endl;
    cout << "rotation mat: "<< endl << dataset.rotations[i_symop];
    cout << "inverse: " << endl << inverse(to_double(dataset.rotations[i_symop]));
    cout << "transition v: " << dataset.translations[i_symop] << endl;
    auto euler = get_Euler_from_sym_matrix_spg(dataset.rotations[i_symop], latt, is_proper);
    printf("Euler angle: %f %f %f, proper? %d\n", euler[0], euler[1], euler[2], is_proper);
    auto Wmat = compute_W_matrix(latt, posi, atom_types, kprime,
                                 dataset.rotations[i_symop],
                                 dataset.translations[i_symop], map_type_abfs, "orig");

    // std::cout << Wmat;
    cout << "Testing unitary property " << endl;
    auto iden = Wmat * transpose(Wmat, true);
    matrix<cplxdb> iden_ref(iden.nr, iden.nc);
    iden_ref.set_diag(1);
    assert(iden == iden_ref);

    map_type_abfs.clear();
}

int main (int argc, char *argv[])
{
    test_prim_nacl_spd();
    test_prim_diamond_spd();
    return 0;
}
