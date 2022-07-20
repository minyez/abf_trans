#include "rotate.h"
#include <cmath>

double * get_rotation_matrix_xyz(double rotmat_spg[3][3])
{
    double * rotmat_xyz = new double [9];
    return rotmat_xyz;
}

std::array<double, 3> get_Euler_from_rotation_matrix(double rotmat[3][3])
{
    std::array<double, 3> euler;
    return euler;
}

void get_Wigner_small_d_matrix_from_Euler(unsigned l, double beta, double* dmat)
{
    auto msize = get_msize(l);
}

void get_Wigner_D_matrix_from_Euler(cplxdb* Dmat, unsigned l, std::array<double, 3> euler_angle)
{
    auto msize = get_msize(l);
}

void get_RSH_Delta_matrix_from_Euler(cplxdb* Delta, unsigned l, std::array<double, 3> euler_angle)
{
    auto msize = get_msize(l);
}
