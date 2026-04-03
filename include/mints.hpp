#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <tuple>
#include "omp.h"

#include "libint2.hpp"
#include "xtensor.hpp"

#include "constants.hpp"



namespace sf {

constexpr double sqrt3_2 = 0.86602540378443859658830206171842291951179504394531;
const xt::xtensor_fixed<double, xt::xshape<3,3>> tf1 = {{{0., 1., 0.},
                                                         {0., 0., 1.},
                                                         {1., 0., 0.}}};

const xt::xtensor_fixed<double, xt::xshape<5,6>> tf2 = {{{     0., 1., 0.,       0., 0., 0.},
                                                         {     0., 0., 0.,       0., 1., 0.},
                                                         {   -0.5, 0., 0.,     -0.5, 0., 1.},
                                                         {     0., 0., 1.,       0., 0., 0.},
                                                         {sqrt3_2, 0., 0., -sqrt3_2, 0., 0.}}};

const xt::xtensor_fixed<double, xt::xshape<7,10>> tf3 = {{{0,                    1.0606601717798213,   0,                    0,                   0,  0,                  -0.79056941504209483,  0,                   0,                  0},
                                                          {0,                    0,                    0,                    0,                   1,  0,                   0,                    0,                   0,                  0},
                                                          {0,                   -0.27386127875258306,  0,                    0,                   0,  0,                  -0.61237243569579452,  0,                   1.0954451150103322, 0},
                                                          {0,                    0,                   -0.67082039324993691,  0,                   0,  0,                   0,                   -0.67082039324993691, 0,                  1},
                                                          {0.61237243569579452,  0,                    0,                   -0.27386127875258306, 0,  1.0954451150103322,  0,                    0,                   0,                  0},
                                                          {0,                    0,                    0.86602540378443865,  0,                   0,  0,                   0,                   -0.86602540378443865, 0,                  0},
                                                          {0.79056941504209483,  0,                    0,                   -1.0606601717798213,  0,  0,                   0,                    0,                   0,                  0}}};

double sfact2(int n);
std::vector<std::array<int, 3>> cart_ordering(int ltot);
std::vector<std::vector<libint2::Shell>> read_g94_basis_library(std::string file_dot_g94,
                                                                bool force_cartesian_d=false,
                                                                bool throw_if_missing=true,
                                                                std::string locale_name=std::string{"POSIX"});
std::vector<libint2::Atom> make_atoms(const std::string& xyzfile);
std::vector<libint2::Shell> make_shells(const std::vector<libint2::Atom>& atoms, bool pure, const std::string& file_dot_g94="../data/gbs/cc-pvdz.g94");
std::vector<libint2::Shell> _make_shells_cart_noembed(const std::vector<libint2::Atom>& atoms, const std::string& file_dot_g94="../data/gbs/cc-pvdz.g94");
std::vector<size_t> get_shell2bf(const std::vector<libint2::Shell>& shells);
size_t get_max_nprim(const std::vector<libint2::Shell>& shells);
int get_lmax(const std::vector<libint2::Shell>& shells);
xt::xtensor<double, 2> get_olp(const std::vector<libint2::Shell>& shells);
xt::xtensor<double, 2> get_kin(const std::vector<libint2::Shell>& shells);
xt::xtensor<double, 2> get_ext(const std::vector<libint2::Shell>& shells, const std::vector<libint2::Atom>& atoms);
std::pair<xt::xtensor<double, 2>, xt::xtensor<double, 2>> get_JK(const std::vector<libint2::Shell>& shells, const xt::xtensor<double, 2>& D);
xt::xtensor<double, 2> get_J(const std::vector<libint2::Shell>& shells, const xt::xtensor<double, 2>& D);
xt::xtensor<double, 2> get_J_par(const std::vector<libint2::Shell>& shells, const xt::xtensor<double, 2>& D);



struct Molecule
{
    Molecule(const std::string& xyzfile, const std::string& name);
    Molecule() = default;
    ~Molecule() = default;
    xt::xtensor<double, 2> make_tf() const;

private:
    double get_e_nuc() const;

public:
    std::vector<libint2::Atom>  atoms;
    std::vector<libint2::Shell> shells_pure;
    std::vector<libint2::Shell> shells_cart;
    std::vector<libint2::Shell> shells_cart_noembed;
    std::string                 xyzfile;
    std::string                 name;
    int                         natom;
    int                         nshell;
    int                         nbf_cart;
    int                         nbf_pure;
    std::vector<int>            nbf_in_shells_cart;
    std::vector<int>            nbf_in_shells_pure;
    std::vector<std::string>    symbols;
    std::vector<double>         rms;
    std::vector<int>            numbers;
    xt::xtensor<double, 2>      xyz;
    double                      e_nuc;
    int                         ne;
    int                         nocc;
};



struct BFs {
    struct BF {
        xt::xtensor<double, 1> exponents;
        xt::xtensor<double, 1> coefficients;
        xt::xtensor<double, 1> normfactors;
        xt::xtensor_fixed<double, xt::xshape<3>> center;
        xt::xtensor_fixed<int, xt::xshape<3>> lxyz;
        int nprim = exponents.size();
    };

    BFs(const Molecule& mol);
    BFs() = default;
    ~BFs() = default;
    xt::xtensor<double, 1> get_ao_val(double x, double y, double z) const;
    xt::xtensor<double, 2> get_ao_val(const xt::xtensor<double, 1>& x, const xt::xtensor<double, 1>& y, const xt::xtensor<double, 1>& z) const;
    std::array<xt::xtensor<double, 1>, 3> get_ao_grad(double x, double y, double z);
    std::array<xt::xtensor<double, 2>, 3> get_ao_grad(const xt::xtensor<double, 1>& x, const xt::xtensor<double, 1>& y, const xt::xtensor<double, 1>& z);
    std::vector<BF> bfs;
    int nbf;
};



} // end namespace sf

