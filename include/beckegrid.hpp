#ifndef INCLUDE_BECKEGRID_HPP
#define INCLUDE_BECKEGRID_HPP

#include <tuple>
#include <vector>
#include <string>
#include <cassert>

#include "xtensor.hpp"

#include "lebedev.hpp"
#include "constants.hpp"


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

double sk(int k, double nuij);


//>! sij to be modified
void sk(int k, const xt::xtensor<double, 1>& nuij, xt::xtensor<double, 1>& sij);


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
 * @class Becke分子网格类
 * @brief 构建Becke网格，网格与分子分离，网格初始化的参数应该由分子类提供
 *
 * @param  symbols   {std::vector<std::string>} - element symbols,       given by `Molecule`
 * @param      rms   {std::vector<double>}      - Bragg-Slater radii,    given by `Molecule`
 * @param numnbers   {std::vector<size_t>}      - atomic number,         given by `Molecule`
 * @param      xyz   {xt::xtensor<double, 2>}   - (N, 3) xyz coordinate, given by `Molecule`
 * @param    natom   {size_t}                   - number of atoms,       given by 'Molecule'
 * 
 * @param nrad  = 75   {size_t}                 - radial grids number, nrad
 * @param nleb  = 29   {size_t}                 - Lebedev quadrature order, THIS IS NOT nang
 * @param k     = 4    {size_t}                 - Becke switching function order
 * @param bised = true {bool}                   - if assign different rms to different elements
 * */
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


    //>! 计算所有原子在某点处的权重
    xt::xtensor<double, 1> get_weight_s(const double x, const double y, const double z);


    /**
     * @brief 计算所有原子在一堆点处的权重，返回 (natom, ngrid)矩阵，每列总和为1
     *        可以用来给平面分子对网格点最大权重画图
    */
    xt::xtensor<double, 2> get_weight_p(const xt::xtensor<double, 2>& grid);


    //>! 计算一个原子在自身原子网格上的权重 grid is atgrid; 内部调用 `get_weight_p` grid shape (natgrid, 3)
    xt::xtensor<double, 1> get_weight(const xt::xtensor<double, 2>& grid, int iatom);


    /**
     * @brief 构建分子网格. 当该函数被调用时未在调用构造函数时初始化的数据成员被初始化
     *        由于 `grid xtensor<double, 4>` 用不到，该函数已被弃用
     *
     * @return grid {xt::xtensor<double, 4>} - (natom, nrad, nang, 3+natom) 前三列是网格点坐标
     * @return grid_global {xt::xtensor<double, 3>} - (natom, natgrid, 3) 只记录坐标
     * */
    [[deprecated("this method returns one unnecessary array, use build_grid2 instead\n")]]
    std::tuple<xt::xtensor<double, 4>, xt::xtensor<double, 3>> build_grid();


    /**
     * @brief 构造分子网格，未在构造函数中初始化的成员变量在该函数内部被初始化
     * @return grid_global {std::vector<xtensor<double, 2>>} - [natom, (natgrid, 3)] 存储各个原子网格的坐标
     */
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

    // 剩下的数据成员在调用构造函数时不被初始化, 当调用 build_grid2 时进行赋值
    xt::xtensor<double, 1>   zcheb; // (nrad, ) 所有原子共用一套
    xt::xtensor<double, 2>   xwleb; // (4, nang) 所有原子共用一套, 前三行是xyz坐标, 最后一行是权重; lebedev on unit sphere
    std::vector<xt::xtensor<double, 2>> xrwcheb; // [natom, (3, nrad)] xrw; chebyshevII radial grid
    std::vector<xt::xtensor<double, 1>> weights; // [natom, natgrid] 各原子在自身网格上的权重 atomic weight wa
};


} // end namespace

#endif
