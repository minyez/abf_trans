#include "../src/abf.h"
#include <cassert>
#include <iostream>

using namespace std;

void test_abf_prim_diamond_tier1()
{
    std::vector<int> atypes;
    atypes.push_back(6);
    atypes.push_back(6);
    map_type_abfs[6].push_back({0, 1, 0});
    map_type_abfs[6].push_back({0, 2, 0});
    map_type_abfs[6].push_back({1, 2, 0});
    map_type_abfs[6].push_back({0, 2, 1});
    map_type_abfs[6].push_back({1, 2, 1});
    map_type_abfs[6].push_back({1, 3, 2});
    ABF basis(atypes, map_type_abfs);
    cout << "Total number of wfc. basis functions: " << basis.get_number_of_total_abfs() << endl;
    assert(basis.get_number_of_total_abfs() == 28);
    map_type_abfs.clear();
}


void test_basbas_prim_diamond_tier1()
{
    std::vector<int> atypes;
    atypes.push_back(6);
    atypes.push_back(6);
    map_type_abfs[6].push_back({0, 1, 0});
    map_type_abfs[6].push_back({0, 2, 0});
    map_type_abfs[6].push_back({0, 3, 0});
    map_type_abfs[6].push_back({0, 4, 0});
    map_type_abfs[6].push_back({0, 5, 0});
    map_type_abfs[6].push_back({0, 6, 0});
    map_type_abfs[6].push_back({0, 7, 0});
    map_type_abfs[6].push_back({0, 8, 0});
    map_type_abfs[6].push_back({0, 1, 1});
    map_type_abfs[6].push_back({0, 2, 1});
    map_type_abfs[6].push_back({0, 3, 1});
    map_type_abfs[6].push_back({0, 4, 1});
    map_type_abfs[6].push_back({0, 5, 1});
    map_type_abfs[6].push_back({0, 6, 1});
    map_type_abfs[6].push_back({0, 2, 2});
    map_type_abfs[6].push_back({0, 3, 2});
    map_type_abfs[6].push_back({0, 4, 2});
    map_type_abfs[6].push_back({0, 5, 2});
    map_type_abfs[6].push_back({0, 6, 2});
    map_type_abfs[6].push_back({0, 7, 2});
    map_type_abfs[6].push_back({0, 3, 3});
    map_type_abfs[6].push_back({0, 4, 3});
    map_type_abfs[6].push_back({0, 5, 3});
    map_type_abfs[6].push_back({0, 4, 4});
    ABF basis(atypes, map_type_abfs);
    cout << "Total number of product basis functions: " << basis.get_number_of_total_abfs() << endl;
    assert(basis.get_number_of_total_abfs() == 172);
    map_type_abfs.clear();
}

int main (int argc, char *argv[])
{
    test_abf_prim_diamond_tier1();
    test_basbas_prim_diamond_tier1();
    return 0;
}
