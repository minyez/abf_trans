#include "rotate.h"
#include "constants.h"
#include "mathtools.h"
#include <iostream>
#include "sph.h"

matrix<double> get_sym_matrix_xyz(const matrix<int> &rotmat_spg,
                                  const matrix<double> &lattice)
{
    assert(rotmat_spg.nc == 3 && rotmat_spg.nr == 3);
    matrix<double> rotmat_xyz(3, 3);
    rotmat_xyz = transpose(lattice) * to_double(rotmat_spg) * inverse(transpose(lattice));
    return rotmat_xyz;
}

matrix<int> sym_matrix_spg_from_Euler(double alpha, double beta, double gamma, bool is_proper,
                                      const matrix<double> &lattice)
{
    matrix<int> rotmat_spg(3, 3);
    rotmat_spg = 0;
    matrix<double> rotmat_xyz = sym_matrix_xyz_from_Euler(alpha, beta, gamma, is_proper);
    auto rotmat_spg_db = inverse(transpose(lattice)) * rotmat_xyz * transpose(lattice);
    for (int i = 0; i < rotmat_spg_db.size(); i++)
        if (rotmat_spg_db.c[i] > 0.9)
        {
            rotmat_spg.c[i] = 1;
        }
        else if (rotmat_spg_db.c[i] < -0.9)
        {
            rotmat_spg.c[i] = 1;
        }
        else
            throw std::invalid_argument("sym_matrix_spg_from_Euler: ");
    if(!is_proper)
        rotmat_spg *= -1;
    return rotmat_spg;
}


matrix<double> sym_matrix_xyz_from_Euler(double alpha, double beta, double gamma, bool is_proper)
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

    if(!is_proper)
        symmat *= -1;

    return symmat;
}

std::array<double, 3> get_Euler_from_sym_matrix_xyz(const matrix<double> &rotmat_xyz,
                                                    bool &is_proper)
{
    assert(rotmat_xyz.nc == 3 && rotmat_xyz.nr == 3);
    double alpha, beta, gamma;
    matrix<double> prop_rot = rotmat_xyz;
    const double zerothres = 1e-8;
    if (rotmat_xyz.det() < 0)
        is_proper = false;
    else
        is_proper = true;
    if (!is_proper)
        prop_rot *= -1;

    // Version 1
    // beta = std::acos(prop_rot(2, 2));
    // double sinbeta = std::sin(beta);
    // /* std::cout << beta << std::endl; */
    // if (std::fabs(sinbeta) > zerothres)
    // {
    //     // handle value slightly larger than 1 or smaller than -1 during division
    //     double divided_by_sinb;
    //     divided_by_sinb = prop_rot(2, 0) / sinbeta;
    //     if (fabs(divided_by_sinb - 1) < zerothres && divided_by_sinb > 0)
    //         divided_by_sinb = 1.0;
    //     if (fabs(divided_by_sinb + 1) < zerothres && divided_by_sinb < 0)
    //         divided_by_sinb = -1.0;
    //     alpha = std::acos(divided_by_sinb);
    //     if (prop_rot(2, 1) < 0)
    //         alpha = 2*PI - alpha;

    //     divided_by_sinb = - prop_rot(0, 2) / sinbeta;
    //     if (fabs(divided_by_sinb - 1) < zerothres && divided_by_sinb > 0)
    //         divided_by_sinb = 1.0;
    //     if (fabs(divided_by_sinb + 1) < zerothres && divided_by_sinb < 0)
    //         divided_by_sinb = -1.0;
    //     gamma = std::acos(divided_by_sinb);
    //     if (prop_rot(1, 2) < 0)
    //         gamma = 2*PI - gamma;
    // }
    // else
    // {
    //     // always set gamma as 0 for this case.
    //     gamma = 0.;
    //     // symmat(1,1) = cos(a+g) if cosb > 0 i.e. symmat(2,2) > 0
    //     //             = cos(a-g) if cosb < 0
    //     alpha = std::acos(prop_rot(1, 1));
    //     if (prop_rot(0, 1) * prop_rot(2, 2) < 0)
    //         alpha = 2*PI - alpha;
    // }

    // Version 2, adapt from aims code
    if (fabs(prop_rot(2, 0)) > zerothres || fabs(prop_rot(2, 1)) > zerothres)
    {
        alpha = std::atan2(prop_rot(2, 1), prop_rot(2, 0));
        if (alpha < 0) alpha += 2.0*PI;
        gamma = std::atan2(prop_rot(1, 2), -prop_rot(0, 2));
        if (gamma < 0) gamma += 2.0*PI;
        if (fabs(prop_rot(2, 0)) > fabs(prop_rot(2, 1)))
            beta = std::atan2(prop_rot(2, 0) / std::cos(alpha), prop_rot(2, 2));
        else
            beta = std::atan2(prop_rot(2, 1) / std::sin(alpha), prop_rot(2, 2));
    }
    else
    {
        // beta =  0: cos(a+g), sin(a+g)
        // beta = pi: cos(pi+a-g), sin(pi+a-g)
        alpha = std::atan2(prop_rot(0, 1), prop_rot(0, 0));
        if (alpha < 0)
            alpha += 2*PI;
        if (prop_rot(2, 2) > 0)
        {
            beta = 0;
            gamma = 0;
        }
        else
        {
            beta = PI;
            gamma = PI;
        }
    }
    return std::array<double, 3>{alpha, beta, gamma};
}

std::array<double, 3> get_Euler_from_sym_matrix_spg(const matrix<int> &rotmat_spg, const matrix<double> &lattice, bool &is_proper)
{
    auto rotmat_xyz = get_sym_matrix_xyz(rotmat_spg, lattice);
    return get_Euler_from_sym_matrix_xyz(rotmat_xyz, is_proper);
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

matrix<cplxdb> get_Wigner_D_matrix_from_Euler(unsigned int l, const std::array<double, 3> &euler_angle, bool is_proper)
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
    if (!is_proper)
        smalld *= std::pow(-1., il);
    return alpha * to_complex(smalld) * gamma;
}

matrix<cplxdb> get_RSH_Delta_matrix_from_Euler(unsigned l, const std::array<double, 3> &euler_angle, bool is_proper, const CODE_CHOICE &rsh_choice)
{
    auto Dmat= get_Wigner_D_matrix_from_Euler(l, euler_angle, is_proper);
    int msize = get_msize(l);
    int il = int(l);
    matrix<cplxdb> Delta(msize, msize);

    for (int m = -il; m <= il; m++)
        for (int mp = -il; mp <= il; mp++)
        {
            if (m == 0 && mp == 0)
                Delta(m+il, mp+il) = Dmat(m+il, mp+il);
            else if (m == 0)
                Delta(m+il, mp+il) = Dmat(il, mp+il) * get_C_matrix_element(mp, mp, rsh_choice) + Dmat(il, -mp+il) * get_C_matrix_element(mp, -mp, rsh_choice);
            else if (mp == 0)
                Delta(m+il, mp+il) = Dmat(m+il, il) * std::conj(get_C_matrix_element(m, m, rsh_choice)) + Dmat(-m+il, il) * std::conj(get_C_matrix_element(m, -m, rsh_choice));
            else
            {
                Delta(m+il, mp+il) = std::conj(get_C_matrix_element(m, m, rsh_choice)) * (Dmat(m+il, mp+il) * get_C_matrix_element(mp, mp, rsh_choice) + Dmat(m+il, -mp+il) * get_C_matrix_element(mp, -mp, rsh_choice)) +
                                     std::conj(get_C_matrix_element(m, -m, rsh_choice)) * (Dmat(-m+il, mp+il) * get_C_matrix_element(mp, mp, rsh_choice) + Dmat(-m+il, -mp+il) * get_C_matrix_element(mp, -mp, rsh_choice));
            }
        }
    // NOTE: filter out small values, mainly due to the last condition above.
    for (int i = 0; i < Delta.size(); i++)
        if (fabs(Delta.c[i]) < 1.e-14)
            Delta.c[i] = 0.0;
    // ensure real matrix
    assert(get_imag(Delta) == 0.0);
    return Delta;
}
