#include "rotate.h"
#include "constants.h"
#include <iostream>

matrix<double> get_sym_matrix_xyz(const matrix<double> &rotmat_spg,
                                       const matrix<double> &lattice)
{
    assert(rotmat_spg.nc == 3 && rotmat_spg.nr == 3);
    matrix<double> rotmat_xyz(3, 3);
    rotmat_xyz = lattice * rotmat_spg * inverse(lattice);
    return rotmat_xyz;
}


matrix<double> sym_matrix_xyz_from_Euler(double alpha, double beta, double gamma, bool add_inversion)
{
    matrix<double> symmat(3, 3);
    const double cosa = std::cos(alpha);
    const double sina = std::sin(alpha);
    const double cosb = std::cos(beta);
    const double sinb = std::sin(beta);
    const double cosg = std::cos(gamma);
    const double sing = std::sin(gamma);

    symmat(0, 0) = cosa*cosb*cosg - sina*sing;
    symmat(0, 1) = sina*cosb*cosg + cosa*sing;
    symmat(0, 2) = - sinb*cosg;
    symmat(1, 0) = - cosa*cosb*sing - sina*cosg;
    symmat(1, 1) = - sina*cosb*sing + cosa*cosg;
    symmat(1, 2) = sinb*sing;
    symmat(2, 0) = cosa*sinb;
    symmat(2, 1) = sina*sinb;
    symmat(2, 2) = cosb;

    if(add_inversion)
        symmat *= -1;

    return symmat;
}

std::array<double, 3> get_Euler_from_sym_matrix_xyz(const matrix<double> &rotmat_xyz,
                                                    bool &is_proper)
{
    assert(rotmat_xyz.nc == 3 && rotmat_xyz.nr == 3);
    double alpha, beta, gamma;
    matrix<double> prop_rot = rotmat_xyz;
    if (rotmat_xyz.det() < 0)
        is_proper = false;
    else
        is_proper = true;
    if (!is_proper)
        prop_rot *= -1;
    beta = std::acos(prop_rot(2, 2));
    double sinbeta = std::sin(beta);
    /* std::cout << beta << std::endl; */
    if (std::fabs(sinbeta) > 1.0e-7)
    {
        // handle value slightly larger than 1 or smaller than -1 during division
        double divided_by_sinb;
        divided_by_sinb = prop_rot(2, 0) / sinbeta;
        if (fabs(divided_by_sinb - 1) < 1e-7 && divided_by_sinb > 0)
            divided_by_sinb = 1.0;
        if (fabs(divided_by_sinb + 1) < 1e-7 && divided_by_sinb < 0)
            divided_by_sinb = -1.0;
        alpha = std::acos(divided_by_sinb);
        if (prop_rot(2, 1) < 0)
            alpha = 2*PI - alpha;

        divided_by_sinb = - prop_rot(0, 2) / sinbeta;
        if (fabs(divided_by_sinb - 1) < 1e-7 && divided_by_sinb > 0)
            divided_by_sinb = 1.0;
        if (fabs(divided_by_sinb + 1) < 1e-7 && divided_by_sinb < 0)
            divided_by_sinb = -1.0;
        gamma = std::acos(divided_by_sinb);
        if (prop_rot(1, 2) < 0)
            gamma = 2*PI - gamma;
    }
    else
    {
        // always set gamma as 0 for this case.
        // symmat(1,1) = cos(a+g) if cosb > 0 i.e. symmat(2,2) > 0
        //             = cos(a-g) if cosb < 0
        alpha = std::acos(prop_rot(1, 1));
        if (prop_rot(0, 1) * prop_rot(2, 2) < 0)
            alpha = 2*PI - alpha;
    }
    return std::array<double, 3>{alpha, beta, gamma};
}
