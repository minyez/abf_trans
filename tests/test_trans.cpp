#include "../src/trans.h"
#include "../src/spglib_utils.h"
#include "testutils.h"

void test_prim_diamond_spd()
{
    using namespace std;
    matrix<double> latt(3, 3), posi(2, 3);
    latt(0, 1) = latt(0, 2) = latt(1, 0) = latt(1, 2) = latt(2, 0) = latt(2, 1) = 0.5;
    posi(1, 0) = posi(1, 1) = posi(1, 2) = 0.25;
    vector<int> atom_types{6, 6};
    SpgDS_c dataset(latt, posi, atom_types);

    map_type_abfs[6].push_back({0});
    map_type_abfs[6].push_back({1});
    map_type_abfs[6].push_back({2});

    vec<double> kprime(3);
    const int i_symop = 0;
    kprime[0] = 0.5;
    cout << "k': "<< kprime << endl;
    cout << "rotation mat: "<< endl << dataset.rotations[i_symop];
    cout << "transition v:" << dataset.translations[i_symop] << endl;
    auto Wmat = compute_W_matrix(latt, posi, atom_types, kprime,
                                 to_double(dataset.rotations[i_symop]),
                                 dataset.translations[i_symop], map_type_abfs);
    std::cout << Wmat;
}

int main (int argc, char *argv[])
{
    test_prim_diamond_spd();
    return 0;
}
