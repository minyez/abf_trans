#include "rotate.h"
#include "constants.h"
#include "mathtools.h"
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

matrix<double> get_Wigner_small_d_matrix_from_Euler_beta(const unsigned int &l, const double &beta)
{
    int msize = get_msize(l);
    int il = int(l);
    double coshb = std::cos(beta/2);
    double sinhb = std::sin(beta/2);
    matrix<double> dmat(msize, msize);
    // TODO: use relations of djm'm to reduce the calculations
    for (int mp = -il; mp <= il; mp++)
        for (int m = -il; m <= il; m++)
        {
            int maxi = std::min(il - mp, il + m);
            int mini = std::max(0, m - mp);
            /* printf("%d %d %d %d\n", mp, m, mini, maxi); */
            for (int i = mini; i <= maxi; i++)
            {
                /* printf("%d %d %d\n", mp, m, i); */
                dmat(mp+il, m+il) += std::pow(-1, i) * std::pow(coshb, 2*il+m-mp-2*i) * std::pow(sinhb, mp-m+2*i)
                    / (factorial(i) * factorial(il-mp-i) * factorial(il+m-i) * factorial(i-m+mp));
            }
            dmat(mp+il, m+il) *= std::pow(-1, mp-m) * std::sqrt(factorial(il+m)*factorial(il-m)*factorial(il+mp)*factorial(il-mp));
        }
    return dmat;
}

matrix<cplxdb> get_Wigner_D_matrix_from_Euler(unsigned int l, const std::array<double, 3> &euler_angle)
{
    int msize = get_msize(l);
    int il = int(l);
    matrix<double> smalld = get_Wigner_small_d_matrix_from_Euler_beta(l, euler_angle[1]);
    matrix<cplxdb> alpha(msize, msize), gamma(msize, msize);
    for (int m = -il; m <= il; m++)
    {
        double angle = m * euler_angle[0];
        alpha(m+il, m+il) = cplxdb{std::cos(angle), -std::sin(angle)};
        angle = m * euler_angle[2];
        gamma(m+il, m+il) = cplxdb{std::cos(angle), -std::sin(angle)};
    }
    return alpha * to_complex(smalld) * gamma;
}
