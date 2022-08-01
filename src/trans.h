#pragma once
#include "abf.h"
#include "cell.h"
#include "rotate.h"

//! compute the transformation matrix W(V;k'). Note kprime is in reciprocal coordiante
matrix<cplxdb> compute_W_matrix(const matrix<double> &lattice, 
                                const matrix<double> &positions, 
                                const vector<int> &atom_types,
                                const vec<double> &kprime,
                                const matrix<double> &rotmat_spg,
                                const vec<double> &transi_spg,
                                const std::map<int, std::vector<abf_id>> &map_type_abfs_in);
