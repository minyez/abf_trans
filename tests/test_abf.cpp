#include "../src/abf.h"
#include <cassert>
#include <iostream>

using namespace std;

void check_search_basis_index(const ABF &basis, int id, int iat, int irf, int l, int m, const CODE_CHOICE &choice)
{
    int temp_iat, temp_irf, temp_l, temp_m;
    basis.get_abf_arlm(id, temp_iat, temp_irf, temp_l, temp_m, choice);

    int temp_id = basis.get_abf_index(temp_iat, temp_irf, temp_m, choice);
    printf("id %d = iat %d irf %d l %d m %d, reverted to id = %d\n", id, temp_iat, temp_irf, temp_l, temp_m, temp_id);

    assert(iat == temp_iat && irf == temp_irf && l == temp_l && m == temp_m );
    assert(temp_id == id);
}

void test_index_arlm_consistency()
{
    const CODE_CHOICE choices[] = {CODE_CHOICE::ORIG, CODE_CHOICE::AIMS, CODE_CHOICE::ABACUS};
    std::vector<int> atypes;
    atypes.push_back(3);
    atypes.push_back(6);
    atypes.push_back(12);
    atypes.push_back(3);
    map_type_abfs[3].push_back({0});
    map_type_abfs[3].push_back({2});
    map_type_abfs[6].push_back({0});
    map_type_abfs[6].push_back({1});
    map_type_abfs[6].push_back({1});
    map_type_abfs[6].push_back({1});
    map_type_abfs[6].push_back({2});
    map_type_abfs[12].push_back({0});
    map_type_abfs[12].push_back({2});
    map_type_abfs[12].push_back({1});
    map_type_abfs[12].push_back({3});
    ABF basis(atypes, map_type_abfs);
    int iat, irf, l, m;
    for (auto choice: choices)
    {
        for (int i = 0; i != basis.get_number_of_total_abfs(); i++)
        {
            basis.get_abf_arlm(i, iat, irf, l, m, choice);
            int i_reverse = basis.get_abf_index(iat, irf, m, choice);
            assert(i == i_reverse);
        }
    }

    map_type_abfs.clear();
}

void test_abf_prim_diamond_tier1()
{
    std::vector<int> atypes;
    atypes.push_back(6);
    atypes.push_back(6);
    map_type_abfs[6].push_back({0});
    map_type_abfs[6].push_back({0});
    map_type_abfs[6].push_back({0});
    map_type_abfs[6].push_back({1});
    map_type_abfs[6].push_back({1});
    map_type_abfs[6].push_back({2});
    ABF basis(atypes, map_type_abfs);
    cout << "Total number of wfc. basis functions: " << basis.get_number_of_total_abfs() << endl;
    /* for (int i = 0; i < atypes.size(); i++) */
    /* { */
    /*     cout << basis.start_index_atom[i] << " " << basis.end_index_atom[i] << endl; */
    /* } */
    assert(basis.get_number_of_total_abfs() == 28);

    check_search_basis_index(basis, 13, 0, 5, 2, 2, CODE_CHOICE::AIMS);
    check_search_basis_index(basis, 14, 1, 0, 0, 0, CODE_CHOICE::AIMS);
    check_search_basis_index(basis, 23, 1, 5, 2, -2, CODE_CHOICE::AIMS);

    map_type_abfs.clear();
}


void test_basbas_prim_diamond_tier1()
{
    std::vector<int> atypes;
    atypes.push_back(6);
    atypes.push_back(6);
    map_type_abfs[6].push_back({0});
    map_type_abfs[6].push_back({0});
    map_type_abfs[6].push_back({0});
    map_type_abfs[6].push_back({0});
    map_type_abfs[6].push_back({0});
    map_type_abfs[6].push_back({0});
    map_type_abfs[6].push_back({0});
    map_type_abfs[6].push_back({0});
    map_type_abfs[6].push_back({1});
    map_type_abfs[6].push_back({1});
    map_type_abfs[6].push_back({1});
    map_type_abfs[6].push_back({1});
    map_type_abfs[6].push_back({1});
    map_type_abfs[6].push_back({1});
    map_type_abfs[6].push_back({2});
    map_type_abfs[6].push_back({2});
    map_type_abfs[6].push_back({2});
    map_type_abfs[6].push_back({2});
    map_type_abfs[6].push_back({2});
    map_type_abfs[6].push_back({2});
    map_type_abfs[6].push_back({3});
    map_type_abfs[6].push_back({3});
    map_type_abfs[6].push_back({3});
    map_type_abfs[6].push_back({4});
    ABF basis(atypes, map_type_abfs);
    cout << "Total number of product basis functions: " << basis.get_number_of_total_abfs() << endl;
    assert(basis.get_number_of_total_abfs() == 172);
    check_search_basis_index(basis, 85, 0, 23, 4, 4, CODE_CHOICE::AIMS);
    check_search_basis_index(basis, 100, 1, 10, 1, -1, CODE_CHOICE::AIMS);
    check_search_basis_index(basis, 122, 1, 16, 2, -2, CODE_CHOICE::AIMS);
    check_search_basis_index(basis, 155, 1, 21, 3, 3, CODE_CHOICE::AIMS);

    map_type_abfs.clear();
}

void test_basbas_prim_lif_tier1()
{
    std::vector<int> atypes;
    atypes.push_back(3);
    atypes.push_back(9);
    map_type_abfs[3].push_back({0});
    map_type_abfs[3].push_back({2});
    map_type_abfs[3].push_back({0});
    map_type_abfs[3].push_back({2});
    map_type_abfs[3].push_back({0});
    map_type_abfs[3].push_back({4});
    map_type_abfs[3].push_back({3});
    map_type_abfs[3].push_back({2});
    map_type_abfs[3].push_back({1});
    map_type_abfs[3].push_back({0});
    map_type_abfs[9].push_back({0});
    map_type_abfs[9].push_back({0});
    map_type_abfs[9].push_back({1});
    map_type_abfs[9].push_back({0});
    map_type_abfs[9].push_back({1});
    map_type_abfs[9].push_back({2});
    map_type_abfs[9].push_back({1});
    map_type_abfs[9].push_back({0});
    
    ABF basis(atypes, map_type_abfs);
    cout << "Total number of product basis functions: " << basis.get_number_of_total_abfs() << endl;
    assert(basis.get_number_of_total_abfs() == 56);

    check_search_basis_index(basis, 36, 0, 8, 1, 1, CODE_CHOICE::AIMS);
    check_search_basis_index(basis, 54, 1, 6, 1, 1, CODE_CHOICE::AIMS);
    check_search_basis_index(basis, 55, 1, 7, 0, 0, CODE_CHOICE::AIMS);

    map_type_abfs.clear();
}

int main (int argc, char *argv[])
{
    test_index_arlm_consistency();
    test_abf_prim_diamond_tier1();
    test_basbas_prim_diamond_tier1();
    test_basbas_prim_lif_tier1();
    return 0;
}
