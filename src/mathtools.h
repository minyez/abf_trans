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
