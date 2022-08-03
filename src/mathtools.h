#include <cmath>

inline double factorial(int n)
{
    if (n > 1)
        return n*factorial(n-1);
    if (n < 0)
        return INFINITY;
    return 1;
}

inline double double_factorial(int n)
{
    if (n > 1)
        return n*double_factorial(n-2);
    if (n < -1)
    {
        if (n%2 == 0)
        {
            return INFINITY;
        }
        else
        {
            return double_factorial(n+2) / (n+2);
        }
    }
    return 1;
}

inline double shift_to_unit(double &v, double lowlim, bool keep_lowlim)
{
    double shift = 0.;
    const double uplim = lowlim + 1.0;
    while (v < lowlim)
    {
        v += 1.0;
        shift -= 1.0;
    }
    while (v > uplim)
    {
        v -= 1.0;
        shift += 1.0;
    }
    if (keep_lowlim && v == uplim)
    {
        v -= 1.0;
        shift += 1.0;
    }
    if (!keep_lowlim && v == lowlim)
    {
        v += 1.0;
        shift -= 1.0;
    }
    return shift;
}
