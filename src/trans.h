#pragma once
#include "abf.h"
#include "cell.h"
#include "rotate.h"

//! compute the transformation matrix M(S;k).
matrix<cplxdb> compute_M_matrix(const matrix<double> &lattice, 
                                const matrix<double> &positions, 
                                const vector<int> &atom_types,
                                const vec<double> &tilde_k,
                                const matrix<int> &rotmat_spg,
                                const vec<double> &transi_spg,
                                const std::map<int, std::vector<abf_id>> &map_type_abfs_in,
                                const CODE_CHOICE & choice);

matrix<cplxdb> compute_M_matrix_aims(const matrix<double> &lattice, 
                                     const matrix<double> &positions, 
                                     const vector<int> &atom_types,
                                     const vec<double> &k,
                                     const matrix<int> &rotmat_spg,
                                     const vec<double> &transi_spg,
                                     const std::map<int, std::vector<abf_id>> &map_type_abfs_in);

matrix<cplxdb> compute_representation_on_equiv_k(const vec<double> &kprime, const matrix<cplxdb> &cmat_at_kprime,
                                                 const matrix<double> &lattice, 
                                                 const matrix<double> &positions, 
                                                 const vector<int> &atom_types,
                                                 const matrix<int> &rotmat_spg,
                                                 const vec<double> &transi_spg,
                                                 const std::map<int, std::vector<abf_id>> &map_type_abfs_in,
                                                 const CODE_CHOICE & choice);
