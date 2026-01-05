#ifndef INCLUDE_BECKEGRID_HPP
#define INCLUDE_BECKEGRID_HPP

#include <tuple>
#include <vector>
#include <string>
#include <cassert>

#include "Eigen/Dense"
#include "xtensor.hpp"

#include "lebedev.hpp"
#include "constants.hpp"

using xxd = Eigen::Matrix<double, -1, -1, Eigen::RowMajor>;
using  xd = Eigen::VectorXd;

namespace beckegrid
{

/// Becke switching functions
inline auto p  = [](const double u) -> double {return 3. / 2 * u - 1. / 2 * u * u * u;};
inline auto f1 = [](const double u) -> double {return p(u);}    ; inline auto s1 = [](const double u) -> double {return 0.5 * (1. - f1(u));};
inline auto f2 = [](const double u) -> double {return p(f1(u));}; inline auto s2 = [](const double u) -> double {return 0.5 * (1. - f2(u));};
inline auto f3 = [](const double u) -> double {return p(f2(u));}; inline auto s3 = [](const double u) -> double {return 0.5 * (1. - f3(u));};
inline auto f4 = [](const double u) -> double {return p(f3(u));}; inline auto s4 = [](const double u) -> double {return 0.5 * (1. - f4(u));};
inline auto f5 = [](const double u) -> double {return p(f4(u));}; inline auto s5 = [](const double u) -> double {return 0.5 * (1. - f5(u));};
inline auto f6 = [](const double u) -> double {return p(f5(u));}; inline auto s6 = [](const double u) -> double {return 0.5 * (1. - f6(u));};


/**
 * @brief 高斯-切比雪夫第二类求积
 * 
 * @param norder=75 {size_t} - 径向节点数目
 * @param rm=1.     {double} - Bragg-Slater 半径
 * 
 * @return zxrw {xt::xtensor<double, 2>} - (4, nrad) z ∈ [1, nrad] x ∈ (1, -1); r ∈ (0, big_radii); w weight (Jacobian included)
*/
xt::xtensor<double, 2> gaussCheby2(const size_t norder=75, const double rm=1.);

/**
 * @brief 有限差分法求解径向泊松方程
 * 
 * @param[in] rho_lm {xt::xtensor<double, 1>} - (nrad, ) 电子密度的球谐展开系数（仅是径向r的函数）
 * @param[in]      l                 {size_t} - 角动量量子数
 * @param[in]  rcheb {xt::xtensor<double, 1>} - (nrad, ) 高斯-切比雪夫第二类求积节点长度
 * @param[in]   rm=1.                {double} - Bragg-Slater radii
 * @param[in]   qn=1.                {double} - partial charge, should be computed via quadrature, NOT atomic number
 * 
 * @return u_lm {xt::xtensor<double, 1>} - (nrad, ) 通过电子密度的球谐展开系数计算原子网格上库伦势的球谐展开系数
*/
xt::xtensor<double, 1> poisson_solver(const xt::xtensor<double, 1>& rho_lm, const size_t l, const xt::xtensor<double, 1>& rcheb, const double rm=1., const double qn=1.);


/**
 * @class Becke分子网格类
 * @brief 构建Becke网格，网格与分子分离，网格初始化的参数应该由分子类提供
 *
 * @param  symbols   {std::vector<std::string>} - 元素符号，应由分子类提供
 * @param      rms   {std::vector<double>}      - Bragg-Slater radii, by `Molecule`
 * @param numnbers   {std::vector<size_t>}      - atomic number, by `Molecule`
 * @param      xyz   {xt::xtensor<double, 2>}   - (N, 3) xyz coordinate, by `Molecule`
 * @param    natom   {size_t}                   - number of atoms, by 'Molecule'
 * @param ncheb=75   {size_t}                   - radial grids number, nrad
 * @param  nleb=29   {size_t}                   - Lebedev quadrature order, THIS IS NOT nang
 * @param      k=4   {size_t}                   - Becke switching function order
 * @param bised=true {bool}                     - if assign different rms to different elements
 * */
class BeckeFuzzyCell
{
public:
    BeckeFuzzyCell() = default;
    BeckeFuzzyCell(const std::vector<std::string>& symbols,
                   const std::vector<double>&          rms,
                   const std::vector<int>&         numbers,
                   const xt::xtensor<double, 2>&       xyz,
                   const size_t                   ncheb=75,
                   const size_t                    nleb=29,
                   const size_t                        k=4,
                   const bool                 biased=true);
    ~BeckeFuzzyCell() = default;

    /// 计算所有原子在某点处的权重
    xt::xtensor<double, 1> get_weight_s(const double x, const double y, const double z);

    /// 计算所有原子在一堆点处的权重，返回 (natom, ngrid)矩阵，每列总和为1
    xt::xtensor<double, 2> get_weight_p(const xt::xtensor<double, 2>& grid);

    /**
     * @brief 构建分子网格. 当该函数被调用时未在调用构造函数时初始化的数据成员被初始化
     *
     * @return grid {xt::xtensor<double, 4>} - (natom, nrad, nang, 3+natom) 前三列是网格点坐标
     * @return grid_global {xt::xtensor<double, 3>} - (natom, natgrid, 3) 只记录坐标
     * */
    std::tuple<xt::xtensor<double, 4>, xt::xtensor<double, 3>> build_grid();
    
    std::vector<std::string> symbols;
    std::vector<double>          rms;
    std::vector<int>         numbers;
    xt::xtensor<double, 2>       xyz;
    size_t                     natom;
    size_t                     ncheb;
    size_t                      nleb;
    size_t                         k;
    bool                      biased;
    // 剩下的数据成员在调用构造函数时不被初始化, 当调用 build_grid 时进行赋值
    xt::xtensor<double, 1>   zcheb; // (ncheb, ) 所有原子共用一套
    xt::xtensor<double, 2>   xcheb; // (natom, ncheb)
    xt::xtensor<double, 2>   rcheb; // (natom, ncheb)
    xt::xtensor<double, 2>   wcheb; // (natom, ncheb); nrad == ncheb
    xt::xtensor<double, 2>   xwleb; // (4, nang) 前三行是xyz坐标, 最后一行是权重
    xt::xtensor<double, 2> weights; // (natom, nradxnang); nleb=29 -> nang=302
};


} // end namespace

#endif