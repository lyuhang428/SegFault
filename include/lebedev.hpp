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
xt::xtensor<double, 2> lebedev_rule(const int n);
}

#endif
