#include "../include/rsh.hpp"

double rsh::rsh(const int l, const int m, const double theta, const double phi) noexcept
{
    assert((l >= 0) && (std::abs(m) <= l));
    const    int        mm = (m < 0 ? -m : m);
    const double     phase = (mm & 1) ? -1. : 1.;
    const double         P = std::sph_legendre(l, mm, theta);
    constexpr double root2 = 1.4142135623730951454746218587388284504414;

    if (m == 0) return P;
    if (m > 0)  return phase * root2 * P * std::cos(mm * phi);
    return phase * root2 * P * std::sin(mm * phi);

    // me
    // const int phase = std::abs(m) % 2 == 0 ? 1 : -1;
    // double res = 0;
    // if (m == 0) res = std::sph_legendre(l, 0, theta);
    // else if (m > 0) res = phase * std::sqrt(2.) * std::sph_legendre(l, m, theta) * std::cos(m * phi);
    // else if (m < 0) res = phase * std::sqrt(2.) * std::sph_legendre(l, -m, theta) * std::sin(-m * phi);
    // return res;
}


