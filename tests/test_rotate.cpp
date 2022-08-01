#include "../src/rotate.h"
#include "../src/constants.h"
#include "testutils.h"
#include <iostream>
using namespace std;

void test_identity()
{
    cout << "Testing identity operation ..." << endl;
    matrix<double> iden(3, 3), latt(3, 3);
    iden.set_diag(1);
    latt.set_diag(5);
    assert(iden == get_sym_matrix_xyz(iden, latt));
    bool is_proper;
    auto euler = get_Euler_from_sym_matrix_xyz(iden, is_proper);
    printf("alpha, beta, gamma = %f %f %f\n", euler[0], euler[1], euler[2]);
    assert(fequal(euler[0], 0.0) && fequal(euler[1], 0.0) && fequal(euler[2], 0.0));
    assert(is_proper);

    cout << get_RSH_Delta_matrix_from_Euler(1, euler, true);
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
    double betas[] = { 0.2, 1.0, PI/3, PI/2, 2*PI/3, 3*PI/4, 5*PI/6 };
    double gammas[] = { 0, PI/4, 5*PI/4, 3*PI/2 };

    for (auto alpha: alphas)
        for (auto beta: betas)
            for (auto gamma: gammas)
            {
                /* printf("testing alpha, beta, gamma = %f %f %f\n", alpha, beta, gamma); */
                symmat = sym_matrix_xyz_from_Euler(alpha, beta, gamma, false);
                /* cout << symmat; */
                euler = get_Euler_from_sym_matrix_xyz(symmat, is_proper);
                /* printf("obtained alpha, beta, gamma = %f %f %f\n", euler[0], euler[1], euler[2]); */
                assert(fequal(euler[0], alpha) && fequal(euler[1], beta) && fequal(euler[2], gamma) && is_proper);

                // with inversion
                symmat = sym_matrix_xyz_from_Euler(alpha, beta, gamma, true);
                /* cout << symmat; */
                euler = get_Euler_from_sym_matrix_xyz(symmat, is_proper);
                /* printf("obtained alpha, beta, gamma = %f %f %f\n", euler[0], euler[1], euler[2]); */
                assert(fequal(euler[0], alpha) && fequal(euler[1], beta) && fequal(euler[2], gamma) && !is_proper);
            }

    /* symmat = sym_matrix_xyz_from_Euler(PI/4, PI/3, PI/6, false); */
    /* cout << symmat; */
    /* euler = get_Euler_from_sym_matrix_xyz(symmat, is_proper); */
    /* printf("alpha, beta, gamma = %f %f %f\n", euler[0], euler[1], euler[2]); */
    /* assert(fequal(euler[0], PI/4) && fequal(euler[1], PI/3) && fequal(euler[2], PI/6) && is_proper); */
}

void test_Wigner_smalld()
{
    double betas[] = {0, 0.2, 1.0, PI/3, PI/2, 2*PI/3, 3*PI/4, 5*PI/6 };
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
                printf("testing RSH Delta for alpha, beta, gamma = %f %f %f\n", alpha, beta, gamma);
                std::array<double, 3> euler{alpha, beta, gamma};
                auto Delta = get_RSH_Delta_matrix_from_Euler(1, euler, true);
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
                Delta = get_RSH_Delta_matrix_from_Euler(1, euler, false);
                assert(Delta == (Delta_ref * cplxdb{-1, 0}));
            }
}

int main (int argc, char *argv[])
{
    test_identity();
    test_inversion();
    test_generate_symmat_nonzero_sinbeta();
    test_Wigner_smalld();
    test_Wigner_D();
    test_RSH_Delta();
    return 0;
}
