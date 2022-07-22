#include "../src/abf.h"
#include <cassert>
#include <iostream>

using namespace std;

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
    assert(basis.get_number_of_total_abfs() == 28);

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

    map_type_abfs.clear();
}

int main (int argc, char *argv[])
{
    test_abf_prim_diamond_tier1();
    test_basbas_prim_diamond_tier1();
    return 0;
}
