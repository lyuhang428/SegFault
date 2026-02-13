#ifndef INCLUDE_DFT_HPP
#define INCLUDE_DFT_HPP

#include <cmath>
#include <iostream>
#include <vector>
#include <chrono>
#include <cassert>

#include "omp.h"
#include "xtensor/containers/xtensor.hpp"
#include "xtensor/views/xview.hpp"
#include "xtensor/generators/xgenerator.hpp"
#include "xtensor/io/xnpy.hpp"
#include "xc/xc.h"
#include "libint2.hpp"

#include "mints.hpp"
#include "beckegrid.hpp"
#include "constants.hpp"
#include "lblas.hpp"


namespace sf::DFT
{
    
    /**
     * @class DFT
     *
     * @param xyzfile              {std::string}               - .xyz 文件
     * @param name                 {std::string}               - .94 基组文件 // .json 格式基组文件
     * @param mol                  {sf::Molecule}              - [类对象] 分子类对象
     * @param becke                {beckegrid::BeckeFuzzyCell} - [类对象] Becke分子网格类对象
     * @param natom                {int}                       - 原子数
     * @param nbf_?                {int}                       - 基函数数量
     * @param ne                   {int}                       - 电子数
     * @param nocc                 {int}                       - ne//2
     * @param nrad                 {int}                       - 径向网格数
     * @param nang                 {int}                       - 角度网格数
     * @param natgrid              {int}                       - 原子网格数 natgrid = nrad * nang
     * @param ngrid                {int}                       - 分子网格数 ngrid = natom * nrad * nang
     * @param kinetic_energies     {std::vector<double>}       - 用来记录每一步迭代能量分解的动能项
     * @param external_energies    {std::vector<double>}       - 用来记录每一步迭代能量分解的外势项
     * @param hartree_energies     {std::vector<double>}       - 用来记录每一步迭代能量分解的库伦项
     * @param exchange_energies    {std::vector<double>}       - 用来记录每一步迭代能量分解的交换项
     * @param correlation_energies {std::vector<double>}       - 用来记录每一步迭代能量分解的相关项
     * @param etots                {std::vector<double>}       - 用来记录每一步迭代总能量
     * @param ne_quads             {std::vector<double>}       - 用来记录每一步迭代通过求积得到的电子数
     * @param aos_vals             {xt::xtensor<double, 3>}    - (nao, natom, natgrid) 基函数在网格上的函数值
     */
    class DFT {
    public:
        //>! 该构造函数只初始化 this->xyzfile, this->basefile. 其余数据成员在调用 DFT::init() 和 DFT::scf() 时初始化
        DFT(const std::string& xyzfile, const std::string& name);
        ~DFT() = default;

        /**
         * @brief 初始化分子网格，分子；计算基函数在网格上的值，// 计算全局距离和固体角；计算实球谐函数；计算单位球面上的球谐基
         *
         * @param radial_points=75      {int}         - 径向网格数，默认75
         * @param angular_level=29      {int}         - 列别杰夫求积阶数，默认29（302）
         * @param k=4                   {int}         - Becke switching function level, default 4
         * @param biased=true           {bool}        - 是否对不同的元素使用不同的 Bragg-Slater radii, default true
         * @param radial_scheme="becke" {std::string} - 使用何种径向网格，默认 becke TODO
         */
        void init(int                radial_points=75,
                  int                angular_level=29,
                  int                k=4,
                  bool               biased=true,
                  const std::string& radial_scheme="becke");

        /**
         * @brief 执行 SCF 迭代计算
         *
         * @param maxiter       = 30     {int}         - 最大迭代次数，默认30
         * @param e_convergence = 1e-8   {double}      - 能量收敛标准，默认1e-6
         * @param d_convergence = 1e-8   {double}      - 密度收敛标准，默认1e-6
         * @param nbuffer       = -1     {int}         - DIIS 残差向量长度，默认-1,使用所有历史残差 TODO
         * @param initial_guess = "core" {std::string} - 密度矩阵初始猜测类型 TODO
         * @param X_id                   {int}         - 交换泛函id, 参考 libxc
         * @param C_id                   {int}         - 相关泛函id, 参考 libxc
         */
        void scf(int                maxiter=30,
                 double             e_convergence=1e-8,
                 double             d_convergence=1e-8,
                 int                nbuffer=-1,
                 const std::string& initial_guess="core",
                 int                X_id=1,
                 int                C_id=8, 
                 bool               pure=true);

    private:
        std::string                      xyzfile;
        std::string                      name;
        sf::Molecule                     mol;   // 含有更多数据成员
        beckegrid::BeckeFuzzyCell        becke; // 含有更多数据成员
        int                              natom;
        int                              nbf_cart;
        int                              nbf_pure;
        int                              ne;
        int                              nocc;
        int                              nrad;
        int                              nang;
        int                              natgrid;
        int                              ngrid;
        std::vector<double>              kinetic_energies;
        std::vector<double>              external_energies;
        std::vector<double>              hartree_energies;
        std::vector<double>              exchange_energies;
        std::vector<double>              correlation_energies;
        std::vector<double>              etots;
        std::vector<double>              ne_quads;
        std::vector<xt::xtensor<double, 2>> aos_vals_cart; // [natom, (nbf_cart, natgrid)]
        std::vector<xt::xtensor<double, 2>> aos_vals_pure; // [natom, (nbf_pure, natgrid)]
        

        //>! 开始 SCF 前输出系统信息
        void header_log() const;


        //>! 输出每一步 SCF 的能量分解信息
        void energy_log() const;


        /**
         * @brief commutative Pulay Direct Inversion of the Iterative Subspace (cDIIS)
         * @param[inout] fock     {xxd}              - 待与旧的 fock 混合
         * @param[in]    focks    {std::vector<xxd>} - 历史 Fock 矩阵
         * @param[in]    diis_res {std::vector<xxd>} - 残差向量
        */
        static void cDIIS(xt::xtensor<double, 2>& fock, const std::vector<xt::xtensor<double, 2>>& focks, const std::vector<xt::xtensor<double, 2>>& diis_res);
    };



} // end namespace

#endif
