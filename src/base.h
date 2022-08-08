#pragma once
#include <complex>
#include <cmath>

enum CODE_CHOICE { ORIG, AIMS };
constexpr const char * CODE_CHOICE_STR[] { "original", "FHI-aims" };
enum KRMODE { R, K };
constexpr const char * KRMODE_STR[] { "R", "K" };

typedef std::complex<double> cplxdb;

constexpr static const double DOUBLE_EQUAL_THRES = 1e-10;

inline size_t get_msize(unsigned l)
{
    return 2*l+1;
}

template <typename T, typename T1, typename T2>
T norm(const T v[], const T1 &n, const T2 power)
{
    T norm = 0;
    for (int i = 0; i < n; i++)
        norm += std::pow(std::fabs(v[i]), power);
    return std::pow(norm, 1./power);
}

template <typename Tv, typename T1, typename T2>
Tv norm(const std::complex<Tv> v[], const T1 &n, const T2 power)
{
    Tv norm = 0;
    for (int i = 0; i < n; i++)
        norm += std::pow(std::fabs(v[i]), power);
    return std::pow(norm, 1./power);
}
