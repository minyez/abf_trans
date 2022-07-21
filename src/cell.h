#pragma once
#include "matrix.h"
#include <vector>
#include <map>
#include <set>

using std::vector;
using std::map;
using std::set;

extern matrix<double> latt;
extern matrix<double> posi_frac;
extern vector<int> types;
extern set<int> inequiv_types;
extern map<int, vector<int>> map_type_iatoms;

void generate_map_type_iatom(const vector<int> &atypes, set<int> &inequiv_types, map<int, vector<int>> &amap);
