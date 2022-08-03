#pragma once
#include <string>
#include <vector>
#include <map>
#include <set>
#include "matrix.h"
#include "abf.h"

using std::string;
using std::vector;
using std::map;
using std::set;

void read_cell(const string &cellfile, matrix<double> &out_latt,
               matrix<double> &out_posi_frac, vector<int> &out_types);

void clean_cell(matrix<double> &in_latt,
                matrix<double> &in_posi_frac, vector<int> &in_types);

void read_abf_ids(const string &abffile, const set<int> &types_verify, map<int, vector<abf_id>> &map_t_abf);

matrix<cplxdb> read_mtx_cplxdb(const string &mtxfile);

matrix<double> read_mtx_double(const string &mtxfile);
