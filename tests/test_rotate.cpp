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

int main (int argc, char *argv[])
{
    test_identity();
    test_inversion();
    test_generate_symmat_nonzero_sinbeta();
    return 0;
}
