#pragma once

#include <tuple>
#include <vector>
#include <string>
#include <cassert>

#include "xtensor.hpp"

#include "lebedev.hpp"
#include "constants.hpp"


namespace beckegrid
{
inline auto p  = [](const double u) -> double {return 3. / 2 * u - 1. / 2 * u * u * u;};
inline auto f1 = [](const double u) -> double {return p(u);}    ; inline auto s1 = [](const double u) -> double {return 0.5 * (1. - f1(u));};
inline auto f2 = [](const double u) -> double {return p(f1(u));}; inline auto s2 = [](const double u) -> double {return 0.5 * (1. - f2(u));};
inline auto f3 = [](const double u) -> double {return p(f2(u));}; inline auto s3 = [](const double u) -> double {return 0.5 * (1. - f3(u));};
inline auto f4 = [](const double u) -> double {return p(f3(u));}; inline auto s4 = [](const double u) -> double {return 0.5 * (1. - f4(u));};
inline auto f5 = [](const double u) -> double {return p(f4(u));}; inline auto s5 = [](const double u) -> double {return 0.5 * (1. - f5(u));};
inline auto f6 = [](const double u) -> double {return p(f5(u));}; inline auto s6 = [](const double u) -> double {return 0.5 * (1. - f6(u));};

double sk(int k, double nuij);
void sk(int k, const xt::xtensor<double, 1>& nuij, xt::xtensor<double, 1>& sij);
xt::xtensor<double, 2> gaussCheby2(const size_t norder=75, const double rm=1.);


class BeckeFuzzyCell
{
public:
    BeckeFuzzyCell() = default;
    BeckeFuzzyCell(const std::vector<std::string>& symbols,
                   const std::vector<double>&          rms,
                   const std::vector<int>&         numbers,
                   const xt::xtensor<double, 2>&       xyz,
                   const size_t                    nrad=75,
                   const size_t                    nleb=29,
                   const size_t                        k=4,
                   const bool                 biased=true);
    ~BeckeFuzzyCell() = default;

    xt::xtensor<double, 1> get_weight_s(const double x, const double y, const double z);
    xt::xtensor<double, 2> get_weight_p(const xt::xtensor<double, 2>& grid);
    xt::xtensor<double, 1> get_weight(const xt::xtensor<double, 2>& grid, int iatom);

    [[deprecated("this method returns one unnecessary array, use build_grid2 instead\n")]]
    std::tuple<xt::xtensor<double, 4>, xt::xtensor<double, 3>> build_grid();
    std::vector<xt::xtensor<double, 2>> build_grid2();

    std::vector<std::string> symbols;
    std::vector<double>          rms;
    std::vector<int>         numbers;
    xt::xtensor<double, 2>       xyz;
    size_t                     natom;
    size_t                      nrad;
    size_t                      nleb;
    size_t                         k;
    bool                      biased;
    xt::xtensor<double, 1>   zcheb; // (nrad, ) 所有原子共用一套
    xt::xtensor<double, 2>   xwleb; // (4, nang) 所有原子共用一套, 前三行是xyz坐标, 最后一行是权重; lebedev on unit sphere
    std::vector<xt::xtensor<double, 2>> xrwcheb; // [natom, (3, nrad)] xrw; chebyshevII radial grid
    std::vector<xt::xtensor<double, 1>> weights; // [natom, natgrid] 各原子在自身网格上的权重 atomic weight wa
};


} // end namespace
