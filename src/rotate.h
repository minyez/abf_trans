#pragma once
#include <array>
#include "base.h"
#include "matrix.h"

matrix<double> get_sym_matrix_xyz(const matrix<double> &rotmat_spg,
                                  const matrix<double> &lattice);

matrix<double> sym_matrix_xyz_from_Euler(double alpha, double beta, double gamma, bool add_inversion);

std::array<double, 3> get_Euler_from_sym_matrix_xyz(const matrix<double> &rotmat_xyz,
                                                    bool &is_proper);

//! A convenient function combining 
std::array<double, 3> get_Euler_from_sym_matrix_spg(const matrix<double> &rotmat_spg,
                                                    const matrix<double> &lattice,
                                                    bool &is_proper);

matrix<double> get_Wigner_small_d_matrix_from_Euler_beta(const unsigned int &l, const double &beta);

matrix<cplxdb> get_Wigner_D_matrix_from_Euler(unsigned l, const std::array<double, 3> &euler_angle, bool is_proper);

matrix<cplxdb> get_RSH_Delta_matrix_from_Euler(unsigned l, const std::array<double, 3> &euler_angle, bool is_proper);
