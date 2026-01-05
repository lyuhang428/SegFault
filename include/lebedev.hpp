#ifndef INCLUDE_LEBEDEV_HPP
#define INCLUDE_LEBEDEV_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>

#include "xtensor.hpp"

#include "constants.hpp"


namespace lebedev
{

/**
 * @brief 计算列别杰夫采样点和对应的权重
 * 
 * @param n    列别杰夫求积阶数 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 35, 41, 47, 53, 59, 65, 71, 77, 83, 89, 95, 101, 107, 113, 119, 125, 131
 * 
 * @return xwleb : xt::xtensor<double, 2>, (4, nanag)
 * 
 *          前三行是xyz坐标，最后一行是权重
*/
xt::xtensor<double, 2> lebedev_rule(const int n);




}

#endif