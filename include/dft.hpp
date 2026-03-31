#pragma once

#include <cmath>
#include <iostream>
#include <vector>
#include <cassert>
#include <format>
#include <tuple>
#ifdef TIMEIT
#include <chrono>
#endif

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
    struct DFT {
        //>! 该构造函数只初始化 this->xyzfile, this->basefile. 其余数据成员在调用 DFT::init() 和 DFT::scf() 时初始化
        DFT(const std::string& xyzfile, const std::string& name);
        ~DFT() = default;

        
        //>! 初始化分子网格，分子；计算基函数及其梯度在网格上的值
        void _init_lda(); // 只计算基函数在网格上的值
        void _init_gga(); // 计算基函数及其梯度在网格上的值
        void _init(int         radial_points=75,    // 径向网格数，默认75
                   int         angular_level=29,    // 列别杰夫求积阶数，默认29（302）
                   int         k=3,                 // Becke switching function level, default 3
                   bool        biased=true,         // 是否对不同的元素使用不同的 Bragg-Slater radii, default true
                   std::string radial_scheme="becke", // 使用何种径向网格. TODO: add Treutler-Ahlrichs grid
                   int xc_family=1); // lda calls _init_lda; gga calls _init_gga


        //>! 执行 SCF 迭代计算
        void scf(int         maxiter=30,
                 double      e_convergence=1e-8,
                 double      d_convergence=1e-8,
                 int         nbuffer=-1, // DIIS 残差向量长度，默认-1,使用所有历史残差 TODO: nbuffer=8
                 std::string initial_guess="core", // TODO: add SAD
                 int         X_id=1,
                 int         C_id=8, 
                 bool        pure=true, 
                 int         radial_points=75,    // 径向网格数，默认75
                 int         angular_level=29,    // 列别杰夫求积阶数，默认29（302）
                 int         k=4,                 // Becke switching function level, default 4
                 bool        biased=true,         // 是否对不同的元素使用不同的 Bragg-Slater radii, default true
                 std::string radial_scheme="becke");

    // private:
        std::string                      xyzfile; // .xyz file
        std::string                      name;    // .gbs
        sf::Molecule                     mol;     // 含有更多数据成员
        beckegrid::BeckeFuzzyCell        becke;   // 含有更多数据成员
        int                              natom;
        int                              nbf_cart;
        int                              nbf_pure;
        int                              ne;
        int                              nocc;                   // ne//2
        int                              nrad;                   // 径向网格数
        int                              nang;                   // 角度网格数
        int                              natgrid;                // 原子网格数 natgrid = nrad * nang
        int                              ngrid;                  // 分子网格数 ngrid = natom * nrad * nang
        std::vector<xt::xtensor<double, 2>> aos_vals_cart;       // [natom, (nbf_cart, natgrid)] 笛卡尔基函数在网格上的函数值
        std::vector<xt::xtensor<double, 2>> aos_vals_pure;       // [natom, (nbf_pure, natgrid)] 纯基函数在网格上的函数值 pure = TF @ cart
        std::vector<xt::xtensor<double, 2>> aos_vals_gradx_cart; // [natom, (nbf_cart, natgrid)] 笛卡尔基函数梯度 x-分量在网格上的值
        std::vector<xt::xtensor<double, 2>> aos_vals_grady_cart; // [natom, (nbf_cart, natgrid)] 笛卡尔基函数梯度 y-分量在网格上的值
        std::vector<xt::xtensor<double, 2>> aos_vals_gradz_cart; // [natom, (nbf_cart, natgrid)] 笛卡尔基函数梯度 z-分量在网格上的值
        std::vector<xt::xtensor<double, 2>> aos_vals_gradx_pure; // [natom, (nbf_pure, natgrid)] 纯基函数梯度 x-分量在网格上的值
        std::vector<xt::xtensor<double, 2>> aos_vals_grady_pure; // [natom, (nbf_pure, natgrid)] 纯基函数梯度 y-分量在网格上的值
        std::vector<xt::xtensor<double, 2>> aos_vals_gradz_pure; // [natom, (nbf_pure, natgrid)] 纯基函数梯度 z-分量在网格上的值
        

        //>! 开始 SCF 前输出系统信息
        inline void header_log() const
        {
            std::cout << std::format("Number of atoms          = {:<d}\n", this->natom);
            std::cout << std::format("Number of electrons      = {:<d}\n", this->ne);
            std::cout << std::format("Number of radial grids   = {:<d}\n", this->nrad);
            std::cout << std::format("Number of angular grids  = {:<d}\n", this->nang);
            std::cout << std::format("Number of total grids    = {:<d}\n", this->ngrid);
            std::cout << std::format("Nuclear repulsion energy = {:<.15f} a.u.\n", this->mol.e_nuc);
            std::cout << std::endl;
        }



        /**
         * @brief commutative Pulay Direct Inversion of the Iterative Subspace (cDIIS)
         * @param[inout] fock     {xtensor<double, 2>}              - 待与旧的 fock 混合
         * @param[in]    focks    {std::vector<xtensor<double, 2>>} - 历史 Fock 矩阵
         * @param[in]    diis_res {std::vector<xtensor<double, 2>>} - 残差向量
        */
        static void cDIIS(xt::xtensor<double, 2>& fock, const std::vector<xt::xtensor<double, 2>>& focks, const std::vector<xt::xtensor<double, 2>>& diis_res);
        
        
        //>! electron density at ONE ATOM
        //>! aos_vals can be pure/cart
        xt::xtensor<double, 1> compute_rho(const xt::xtensor<double, 2>& Puv, const xt::xtensor<double, 2>& aos_vals);


        //>! eden and eden grad at ONE ATOM
        //>! aos_vals can be pure/cart
        std::array<xt::xtensor<double, 1>, 3> compute_rho_grad(const xt::xtensor<double, 2>& Puv, 
                                                               const xt::xtensor<double, 2>& aos_vals, 
                                                               const xt::xtensor<double, 2>& aos_vals_gradx, 
                                                               const xt::xtensor<double, 2>& aos_vals_grady, 
                                                               const xt::xtensor<double, 2>& aos_vals_gradz);


        /**
        * @brief numerical quadrature compute XC matrix element of ONE ATOM
        * @remark the first term of gga XC matrix element is also this
        * @param[in] mweighted_prop {xtensor<double, 1>} - prop can be `vx`, `vc`, `ex`, `ec` with shape (natom, natgrid)
        *                                                  `mweights` has shape (natom, natgrid)
        *                                                  mweighted_prop = mweights[iatom] * prop[iatom] ; shape (natgrid, )
        * @param[in] aos_vals {xtensor<double, 2>} - shape (nbf, natgrid), can be pure or cart
        * @return    mat_local {xtensor<double, 2>} - shape (nbf, nbf) <Φ|V|Φ>
        */
        xt::xtensor<double, 2> xc_quadrature_lda(const xt::xtensor<double, 1>& mweighted_prop, const xt::xtensor<double, 2>& aos_vals);
        
        //>! Exuv and Ecuv are computed via xc_quadrature_lda
        //>! this method adds extra term to Kuv and Cuv
        xt::xtensor<double, 2> xc_quadrature_gga(const xt::xtensor<double, 1>& mweighted_vsigma, const xt::xtensor<double, 2>& aos_vals, 
                                                 const xt::xtensor<double, 2>& aos_vals_gradx, 
                                                 const xt::xtensor<double, 2>& aos_vals_grady, 
                                                 const xt::xtensor<double, 2>& aos_vals_gradz, 
                                                 const xt::xtensor<double, 1>& rho_gradx, 
                                                 const xt::xtensor<double, 1>& rho_grady, 
                                                 const xt::xtensor<double, 1>& rho_gradz);


        //>! Puv is density matrix ; op(erator) can be tij, vij, Juv, Exuv, Ecuv, etc.
        double energy_decomposition(const xt::xtensor<double, 2>& Puv, const xt::xtensor<double, 2>& op);
    }; // end struct DFT


} // end namespace

