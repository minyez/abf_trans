#pragma once
#include <array>
#include "base.h"

double * get_rotation_matrix_xyz(double rotmat_spg[3][3]);

std::array<double, 3> get_Euler_from_rotation_matrix(double rotmat[3][3]);

void get_Wigner_small_d_matrix_from_Euler(unsigned l, double beta, double* dmat);

void get_Wigner_D_matrix_from_Euler(cplxdb* Dmat, unsigned l, std::array<double, 3> euler_angle);

void get_RSH_Delta_matrix_from_Euler(cplxdb* Delta, unsigned l, std::array<double, 3> euler_angle);
