#include "../src/rotate.h"
#include "../src/constants.h"
#include "testutils.h"
#include <iostream>
using namespace std;

void test_corotate_k_R()
{
    // the inner product of recip space k and real space R should be invariant
    // when rotating k and R simultaneously
    matrix<int> rot(3, 3);
    rot = 0;
    rot(0, 0) = 1;
    rot(1, 2) = -1;
    rot(2, 1) = 1;
    matrix<double> latt(3, 3);
    latt = 0;
    // a non-standard lattice;
    latt(0, 1) = latt(0, 2) = latt(1, 0) = latt(1, 2) = 1.0;
    latt(2, 0) = latt(2, 1) = 1.0;
    auto recip_latt = transpose(inverse(latt)); // wo 2PI, only as a scaling factor
    cout << latt;
    cout << "Determinant of lattice = " << latt.det() << endl;
    cout << "Determinant of rot.mat = " << to_double(rot).det() << endl;
    auto AAT = latt * transpose(latt);
    auto invAAT = inverse(AAT);
    vec<double> k(3), R(3);
    k[0] = 0.5, k[1] = 0.2, k[2] = -0.1;
    R[0] = 2.5, R[1] = -1.5, R[2] = -0.8;

    double prod = dot(k, R);
    auto R_rot = to_double(rot) * R;
    vec<double> k_rot(3);

    // k_rot = transpose(recip_latt) * k;
    // k_rot = transpose(latt) * to_double(rot) * inverse(transpose(latt)) * k_rot;
    // k_rot = inverse(transpose(recip_latt)) * k_rot;
    k_rot = inverse(transpose(to_double(rot))) * k;
    // cout << AAT;
    // cout << rot;
    // cout << invAAT;
    // cout << (AAT * to_double(rot) * invAAT);
    // cout << (AAT * to_double(rot) * invAAT).det() << endl;
    auto prod_rot = dot(k_rot, R_rot);

    cout << dot(R, R) << " " << dot(R_rot, R_rot) << endl;
    cout << dot(k, k) << " " << dot(k_rot, k_rot) << endl;
    cout << "R(" << R << ") . k(" << k << ") = " << prod << endl;
    cout << dot(transpose(latt) * R, transpose(recip_latt) * k) << endl;
    cout << "R'(" << R_rot << ") . k'(" << k_rot << ") = " << prod_rot << endl;
}

void test_rotate_k()
{
    matrix<double> latt(3, 3);
    latt = 0;
    // a non-standard lattice;
    latt(0, 1) = latt(0, 2) = latt(1, 0) = latt(1, 2) = 1.0;
    latt(2, 0) = latt(2, 1) = 1.0;

    matrix<int> rotmat_spg(3, 3);
    rotmat_spg(0, 1) = rotmat_spg(1, 2) = rotmat_spg(2, 0) = 1;
    vec<double> k(3);
    k[0] = 0.5, k[1] = 0.2, k[2] = -0.1;
    const auto recip_latt = transpose(inverse(latt));
    const auto k_xyz = transpose(recip_latt) * k;
    const auto Vk_xyz = transpose(recip_latt) * rotate_k(rotmat_spg, k, latt);
    const auto rotmat_xyz = transpose(latt) * to_double(rotmat_spg) * inverse(transpose(latt));
    const auto Vk_xyz_direct= rotmat_xyz * k_xyz;
    cout << Vk_xyz << endl;
    cout << Vk_xyz_direct << endl;
    assert(Vk_xyz == Vk_xyz_direct);
}

void test_identity()
{
    cout << "Testing identity operation ..." << endl;
    matrix<int> iden(3, 3);
    matrix<double> latt(3, 3);
    iden.set_diag(1);
    latt.set_diag(5);
    assert(to_double(iden) == get_sym_matrix_xyz(iden, latt));
    bool is_proper;
    auto euler = get_Euler_from_sym_matrix_xyz(to_double(iden), is_proper);
    printf("alpha, beta, gamma = %f %f %f\n", euler[0], euler[1], euler[2]);
    assert(fequal(euler[0], 0.0) && fequal(euler[1], 0.0) && fequal(euler[2], 0.0));
    assert(is_proper);

    // cout << get_RSH_Delta_matrix_from_Euler(1, euler, true, CODE_CHOICE::ORIG);
}

void test_inversion()
{
    cout << "Testing inversion operation ..." << endl;
    matrix<double> inv(3, 3), latt(3, 3);
    inv.set_diag(-1);
    bool is_proper;
    auto euler = get_Euler_from_sym_matrix_xyz(inv, is_proper);
    printf("alpha, beta, gamma = %f %f %f\n", euler[0], euler[1], euler[2]);
    assert(fequal(euler[0], 0.0) && fequal(euler[1], 0.0) && fequal(euler[2], 0.0));
    assert(!is_proper);
}

void test_generate_symmat_nonzero_sinbeta()
{
    double alpha, beta, gamma;
    bool is_proper;
    matrix<double> symmat;
    array<double, 3> euler;

    double alphas[] = { 0, PI/4, 5*PI/4, 3*PI/2 };
    double betas[] = { 0.2, 1.0, PI/3, PI/2, 2*PI/3, 3*PI/4, 5*PI/6};
    double gammas[] = { 0, PI/4, 5*PI/4, 3*PI/2 };

    for (auto alpha: alphas)
        for (auto beta: betas)
            for (auto gamma: gammas)
            {
                /* printf("testing alpha, beta, gamma = %f %f %f\n", alpha, beta, gamma); */
                symmat = sym_matrix_xyz_from_Euler(alpha, beta, gamma, true);
                /* cout << symmat; */
                euler = get_Euler_from_sym_matrix_xyz(symmat, is_proper);
                /* printf("obtained alpha, beta, gamma = %f %f %f\n", euler[0], euler[1], euler[2]); */
                assert(fequal(euler[0], alpha) && fequal(euler[1], beta) && fequal(euler[2], gamma) && is_proper);

                // with inversion
                symmat = sym_matrix_xyz_from_Euler(alpha, beta, gamma, false);
                /* cout << symmat; */
                euler = get_Euler_from_sym_matrix_xyz(symmat, is_proper);
                /* printf("obtained alpha, beta, gamma = %f %f %f\n", euler[0], euler[1], euler[2]); */
                assert(fequal(euler[0], alpha) && fequal(euler[1], beta) && fequal(euler[2], gamma) && !is_proper);
            }

    symmat = sym_matrix_xyz_from_Euler(0, 0, 0, true);
    euler = get_Euler_from_sym_matrix_xyz(symmat, is_proper);
    printf("obtained alpha, beta, gamma = %f %f %f\n", euler[0], euler[1], euler[2]);
    assert(fequal(euler[0], 0.) && fequal(euler[1], 0.) && fequal(euler[2], 0.) && is_proper);

    symmat = sym_matrix_xyz_from_Euler(PI, PI, PI, true);
    euler = get_Euler_from_sym_matrix_xyz(symmat, is_proper);
    printf("obtained alpha, beta, gamma = %f %f %f\n", euler[0], euler[1], euler[2]);
    assert(fequal(euler[0], PI) && fequal(euler[1], PI) && fequal(euler[2], PI) && is_proper);
    /* symmat = sym_matrix_xyz_from_Euler(PI/4, PI/3, PI/6, false); */
    /* cout << symmat; */
    /* euler = get_Euler_from_sym_matrix_xyz(symmat, is_proper); */
    /* printf("alpha, beta, gamma = %f %f %f\n", euler[0], euler[1], euler[2]); */
    /* assert(fequal(euler[0], PI/4) && fequal(euler[1], PI/3) && fequal(euler[2], PI/6) && is_proper); */
}

void test_Wigner_smalld()
{
    double betas[] = {0, 0.2, 1.0, PI/3, PI/2, 2*PI/3, 3*PI/4, 5*PI/6, PI };
    // verify l = 1, the result in p67 and p72 of RoseME57
    // Note that the indices therein are 1, 0, -1, while -1, 0, 1 are used here,
    // following FHI-aims
    for (auto beta: betas)
    {
        matrix<double> smalld1(3, 3);
        double cosb = std::cos(beta);
        double sinb = std::sin(beta);
        smalld1(0, 0) = 0.5*(1+cosb);
        smalld1(0, 1) = sinb / std::sqrt(2);
        smalld1(0, 2) = 0.5*(1-cosb);
        smalld1(1, 0) = - smalld1(0, 1);
        smalld1(1, 1) = cosb;
        smalld1(1, 2) = smalld1(0, 1);
        smalld1(2, 0) = smalld1(0, 2);
        smalld1(2, 1) = smalld1(1, 0);
        smalld1(2, 2) = smalld1(0, 0);
        cout << "Check d(l=1) for beta: " << beta << endl;
        auto smalld1_to_verify = get_Wigner_small_d_matrix_from_Euler_beta(1, beta);
        /* cout << smalld1 << smalld1_to_verify; */
        assert(smalld1 == smalld1_to_verify);
    }
}

void test_Wigner_D()
{

}

void test_RSH_Delta()
{
    // verify l = 1, Eq 53 of BlancoMA97
    double alphas[] = { 0, PI/4, 5*PI/4, 3*PI/2 };
    double betas[] = { 0, 0.2, 1.0, PI/3, PI/2, 2*PI/3, 3*PI/4, 5*PI/6 };
    double gammas[] = { 0, PI/4, 5*PI/4, 3*PI/2 };
    for (auto alpha: alphas)
        for (auto beta: betas)
            for (auto gamma: gammas)
            {
                // printf("testing RSH Delta for alpha, beta, gamma = %f %f %f\n", alpha, beta, gamma);
                std::array<double, 3> euler{alpha, beta, gamma};
                auto Delta = get_RSH_Delta_matrix_from_Euler(1, euler, true, CODE_CHOICE::ORIG);
                matrix<cplxdb> Delta_ref(3, 3);
                const double ca = std::cos(alpha);
                const double cb = std::cos(beta);
                const double cg = std::cos(gamma);
                const double sa = std::sin(alpha);
                const double sb = std::sin(beta);
                const double sg = std::sin(gamma);
                Delta_ref(0, 0) = ca*cg - sa*sg*cb;
                Delta_ref(0, 1) = sa*sb;
                Delta_ref(0, 2) = ca*sg + sa*cg*cb;
                Delta_ref(1, 0) = sg*sb;
                Delta_ref(1, 1) = cb;
                Delta_ref(1, 2) = -cg*sb;
                Delta_ref(2, 0) = -ca*sg*cb - sa*cg;
                Delta_ref(2, 1) = ca*sb;
                Delta_ref(2, 2) = ca*cg*cb - sa*sg;
                /* cout << Delta << Delta_ref; */
                assert(Delta == Delta_ref);
                Delta = get_RSH_Delta_matrix_from_Euler(1, euler, false, CODE_CHOICE::ORIG);
                assert(Delta == (Delta_ref * cplxdb{-1, 0}));
            }
}

void test_Delta_matrix_aims()
{
    matrix<cplxdb> cmat;
    std::array<double, 3> euler;
    const double sqrt3 = sqrt(3);
    const double sqrt5 = sqrt(5);
    const double sqrt6 = sqrt(6);
    const double sqrt10 = sqrt(10);
    {
        euler = { 2*PI, PI/2, 0};
        cplxdb matval_l1[] = { 1, 0, 0,
                               0, 0, 1,
                               0, -1, 0 };
        matrix<cplxdb> mat_l1(3, 3, matval_l1);
        assert(get_RSH_Delta_matrix_from_Euler(1, euler, true, CODE_CHOICE::AIMS) == mat_l1);

        cplxdb matval_l2[] = { 0, 1, 0, 0, 0,
                              -1, 0, 0, 0, 0,
                               0, 0, -0.5, 0, sqrt3*0.5,
                               0, 0, 0, -1, 0,
                               0, 0, sqrt3*0.5, 0, 0.5 };
        matrix<cplxdb> mat_l2(5, 5, matval_l2);
        assert(get_RSH_Delta_matrix_from_Euler(2, euler, true, CODE_CHOICE::AIMS) == mat_l2);

        cplxdb matval_l3[] = { 0.25, 0, 0.25*sqrt3*sqrt5, 0, 0, 0, 0,
                               0, -1, 0, 0, 0, 0, 0,
                            0.25*sqrt3*sqrt5, 0, -0.25, 0, 0, 0, 0,
                               0, 0, 0, 0, -0.25*sqrt6, 0, 0.25*sqrt10,
                               0, 0, 0, 0.25*sqrt6, 0, -0.25*sqrt10, 0,
                               0, 0, 0, 0, 0.25*sqrt10, 0, 0.25*sqrt6,
                               0, 0, 0, -0.25*sqrt10, 0, -0.25*sqrt6, 0};
        matrix<cplxdb> mat_l3_ref(7, 7, matval_l3);
        cout << mat_l3_ref << endl;
        auto mat_l3 = get_RSH_Delta_matrix_from_Euler(3, euler, true, CODE_CHOICE::AIMS);
        cout << mat_l3;
        assert(mat_l3 == mat_l3_ref);
    }
}


int main (int argc, char *argv[])
{
    test_corotate_k_R();
    test_rotate_k();
    test_identity();
    test_inversion();
    test_generate_symmat_nonzero_sinbeta();
    test_Wigner_smalld();
    test_Wigner_D();
    test_RSH_Delta();
    test_Delta_matrix_aims();
    return 0;
}
