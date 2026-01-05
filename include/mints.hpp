#ifndef INCLUDE_MINTS_HPP
#define INCLUDE_MINTS_HPP

#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <tuple>

#include "libint2.hpp"
#include "Eigen/Dense"
#include "xtensor.hpp"
// #include "xtensor/io/xnpy.hpp"

#include "constants.hpp"

using xxd = Eigen::Matrix<double, -1, -1, Eigen::RowMajor>;
using  xd = Eigen::VectorXd;


namespace sf {

constexpr double sqrt3_2 = 0.86602540378443859658830206171842291951179504394531;

//>! cartesian to pure transform matrix
const xxd tf1 = xxd{{{0., 1., 0.}, // -1
                     {0., 0., 1.}, // 0
                     {1., 0., 0.}}}; // 1

const xxd tf2 = xxd{{{     0., 1., 0.,       0., 0., 0.}, // -2
                     {     0., 0., 0.,       0., 1., 0.}, // -1
                     {   -0.5, 0., 0.,     -0.5, 0., 1.}, // 0
                     {     0., 0., 1.,       0., 0., 0.}, // 1
                     {sqrt3_2, 0., 0., -sqrt3_2, 0., 0.}}}; // 2

const xxd tf3 = xxd{{{0,                    1.0606601717798213,   0,                    0,                   0,  0,                  -0.79056941504209483,  0,                   0,                  0}, // -3
                     {0,                    0,                    0,                    0,                   1,  0,                   0,                    0,                   0,                  0}, // -2
                     {0,                   -0.27386127875258306,  0,                    0,                   0,  0,                  -0.61237243569579452,  0,                   1.0954451150103322, 0}, // -1
                     {0,                    0,                   -0.67082039324993691,  0,                   0,  0,                   0,                   -0.67082039324993691, 0,                  1}, //  0
                     {0.61237243569579452,  0,                    0,                   -0.27386127875258306, 0,  1.0954451150103322,  0,                    0,                   0,                  0}, //  1
                     {0,                    0,                    0.86602540378443865,  0,                   0,  0,                   0,                   -0.86602540378443865, 0,                  0}, //  2
                     {0.79056941504209483,  0,                    0,                   -1.0606601717798213,  0,  0,                   0,                    0,                   0,                  0}}}; // 3

double sfact2(int n);

std::vector<std::array<int, 3>> cart_ordering(int ltot);


/**
* @brief reads in all basis sets from a Gaussian94-formatted basis set file (see https://bse.pnl.gov/bse/portal)
*        修改为不将归一化常数合并入缩并系数，后续手动压进时再更改为 true
*        copied from libint2::Basis.h
*
* @param[in] file_dot_g94            {std::string} - file name
* @param[in] force_cartesian_d=false {bool}        - force use of Cartesian d shells, if true
* @param[in] locale_name="POSIX"     {std::string} - specifies the locale to use
*
* @details `export LIBINT_DATA_PATH=/usr/local/share/libint/2.11.2/basis/` if use libint2::Basis{...}
*
* @return vector of basis sets for each element
*/
std::vector<std::vector<libint2::Shell>> read_g94_basis_library(std::string file_dot_g94,
                                                                bool force_cartesian_d=false,
                                                                bool throw_if_missing=true,
                                                                std::string locale_name=std::string{"POSIX"});


std::vector<libint2::Atom> make_atoms(const std::string& xyzfile);


std::vector<libint2::Shell> make_shells(const std::vector<libint2::Atom>& atoms, bool pure, const std::string& file_dot_g94="../data/gbs/cc-pvdz.g94");

//>! cart and DO NOT embed normalization factor into coefficient. This function is only used in extracting g94 info, DO NOT USE IT IN CALC
//>! libint2::do_enforce_unit_normalization 对该函数没有影响，因为 embed=false, 无论什么样的归一化常数都不会乘进系数中，基组信息 will be returned as it is
std::vector<libint2::Shell> _make_shells_cart_noembed(const std::vector<libint2::Atom>& atoms, const std::string& file_dot_g94="../data/gbs/cc-pvdz.g94");

//>! 从 shell 的索引得到 bf 的索引
std::vector<size_t> get_shell2bf(const std::vector<libint2::Shell>& shells);


//>! 从 shell 向量中提取最大缩并数 ； i.e. shells[0] has 8 primitive fn, shells[-1] has 3, ... Then return 8
size_t get_max_nprim(const std::vector<libint2::Shell>& shells);


int get_lmax(const std::vector<libint2::Shell>& shells);


xt::xtensor<double, 2> get_olp(const std::vector<libint2::Shell>& shells);


xt::xtensor<double, 2> get_kin(const std::vector<libint2::Shell>& shells);


xt::xtensor<double, 2> get_ext(const std::vector<libint2::Shell>& shells, const std::vector<libint2::Atom>& atoms);


//>! no optimization, dump TOTAL ERI, just for testing. Production run should not store ERI in mem, do it on the fly
// xt::xtensor<double, 4> get_eri(const std::vector<libint2::Shell>& shells);

/**
 * @brief 不存储整个 ERI ; 在每一步 SCF 动态计算 N^2 的库伦和交换矩阵
 *        copied from github.com/libint/test/
 * @param[in] D {xxd} - 密度矩阵
 * @return {std::pair<xtensor<2>, xtensor<2>>} - 返回库伦 J 和交换 K 矩阵
 * @todo 任务级并行 libint2::Engine 不是线程安全的，每个线程需要创建单独的积分引擎
*/
std::pair<xt::xtensor<double, 2>, xt::xtensor<double, 2>> get_JK(const std::vector<libint2::Shell>& shells, const xxd& D);



struct Molecule
{
    //>! "mol.xyz" ; "../data/gbs/basis.g94"
    Molecule(const std::string& xyzfile, const std::string& name);
    Molecule() = default;
    ~Molecule() = default;

    //>! build cartesian to pure matrix
    xxd make_tf() const;
    xt::xtensor<double, 2> make_tf_xt() const;

private:
    double get_e_nuc() const;

public:
    std::vector<libint2::Atom>        atoms;
    std::vector<libint2::Shell> shells_pure;
    std::vector<libint2::Shell> shells_cart;
    std::vector<libint2::Shell> shells_cart_noembed; // this is only used in reading in g94, DO NOT USE IT IN CALC
    std::string                     xyzfile;
    std::string                        name;
    int                               natom;
    int                              nshell; // same for both cart and pure
    int                            nbf_cart; // total nbf in cart
    int                            nbf_pure; // total nbf in pure, different from nbf_cart if >= d orbitals involved
    std::vector<int>     nbf_in_shells_cart; // i.e. SSSPPD => {1 1 1 3 3 6}
    std::vector<int>     nbf_in_shells_pure; // SSSPPD => {1 1 1 3 3 5} this vector is needed to build cart2pure transform matrix
    std::vector<std::string>        symbols;
    std::vector<double>                 rms; // Bragg-Slater radii
    std::vector<int>                numbers; // charges
    xt::xtensor<double, 2>              xyz; // taken from this->atoms, in Bohr
    double                            e_nuc; // return from this->get_e_nuc()
    int                                  ne; // number of electron
    int                                nocc; // occupation (=ne//2)
};



//>! 用于计算基函数在网格上的函数值
/**
 * @class BFs
 * 
 * @brief All shells of one mol
 * Each shell may contain multiple basis functions (e.g. D has 6 for cartesian, 5 for pure)
 * Each basis function has its own exponent, coefficient, normalization factor, center, and (lx, ly, lz)
 * This class will be used to compute CARTESIAN basis function values on grids
 * And later be converted to SPHERICAL basis function values on grids via conversion matrices (linear combination)
 * 
 * @param nbf {int} - CARTESIAN basis function number
*/
struct BFs {
    //>! BF really plays the game, BFs is merely a collection of BF
    // loop all shells in a mol, then loop all (l,m) combinations in cartesian order to make BF
    struct BF {
        xt::xtensor<double, 1> exponents;
        xt::xtensor<double, 1> coefficients;
        xt::xtensor<double, 1> normfactors;
        xt::xtensor_fixed<double, xt::xshape<3>> center;
        xt::xtensor_fixed<int, xt::xshape<3>> lxyz; // (lx, ly, lz)
        int nprim = exponents.size();
    };

    BFs(const Molecule& mol);
    BFs() = default;
    ~BFs() = default;
    
    //>! all bf value at one point {nbf, } cart
    xt::xtensor<double, 1> ao_val(double x, double y, double z) const;

    //>! all bf values at multiple points {nbf, npoints} cart
    xt::xtensor<double, 2> ao_val(const xt::xtensor<double, 1>& x, const xt::xtensor<double, 1>& y, const xt::xtensor<double, 1>& z) const;

    std::vector<BF> bfs;
    int nbf; // 一定是 cartesian
};



} // end namespace sf

#endif
