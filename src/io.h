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

extern CODE_CHOICE code_choice;
extern KRMODE krmode;

CODE_CHOICE parse_code_choice(const string &cc_in);
KRMODE parse_krmode(const string &mode_in);

double decode_fraction(const string& frac_str);

void read_cell(const string &cellfile, matrix<double> &out_latt,
               matrix<double> &out_posi_frac, vector<int> &out_types);

void clean_cell(matrix<double> &in_latt,
                matrix<double> &in_posi_frac, vector<int> &in_types);

void read_abf_ids(const string &abffile, const set<int> &types_verify, map<int, vector<abf_id>> &map_t_abf);

matrix<cplxdb> read_cplxdb(const string &cmatfile, const string &format_in = "none");
matrix<cplxdb> read_mtx_cplxdb(const string &mtxfile);
matrix<cplxdb> read_csc_cplxdb(const string &cscfile);

matrix<double> read_mtx_double(const string &mtxfile);

void write_mtx_cplxdb(const matrix<cplxdb> &mat, const string &mtxfile,
                      const string &comment = "",
                      double threshold = 1.e-15, bool row_first = true);

void read_matrix_inputs(const string &mat_inputs_fn, std::array<int, 3> &ngs,
                        vector<vec<double>> &krpoints,
                        vector<string> &cmatfns,
                        vector<matrix<cplxdb>> &matrices);

void clear_matrix_inputs(vector<vec<double>> &krpoints, vector<string> &cmatfns,
                         vector<matrix<cplxdb>> &matrices);

