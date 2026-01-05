#ifndef INCLUDE_RSH_HPP
#define INCLUDE_RSH_HPP

#include <cmath>
#include <cassert>

#include "constants.hpp"

namespace rsh
{

/**
 * @brief 计算实球谐函数（Condon-Shortley phase included）
 * 
 * @param l : int, l >= 0 角动量子数
 * @param m : int -l <= m <= l
 * @param theta, phi {double} - solid angle in radian
 * 
 * @return 实球谐函数
*/
double rsh(int l, int m, double theta, double phi) noexcept;

}

#endif