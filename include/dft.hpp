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
        DFT(const std::string& xyzfile, const std::string& name);
        ~DFT() = default;

        void _init_lda();
        void _init_gga();
        void _init(int         radial_points=75,
                   int         angular_level=29,
                   int         k=3,
                   bool        biased=true,
                   std::string radial_scheme="becke",
                   int xc_family=1);

        void scf(int         maxiter=30,
                 double      e_convergence=1e-8,
                 double      d_convergence=1e-8,
                 int         nbuffer=-1,
                 std::string initial_guess="core",
                 int         X_id=1,
                 int         C_id=8, 
                 bool        pure=true, 
                 int         radial_points=75,
                 int         angular_level=29,
                 int         k=3,
                 bool        biased=true,
                 std::string radial_scheme="becke");

        std::string                      xyzfile;
        std::string                      name;
        sf::Molecule                     mol;
        beckegrid::BeckeFuzzyCell        becke;
        int                              natom;
        int                              nbf_cart;
        int                              nbf_pure;
        int                              ne;
        int                              nocc;
        int                              nrad;
        int                              nang;
        int                              natgrid;
        int                              ngrid;
        std::vector<xt::xtensor<double, 2>> aos_vals_cart;
        std::vector<xt::xtensor<double, 2>> aos_vals_pure;
        std::vector<xt::xtensor<double, 2>> aos_vals_gradx_cart;
        std::vector<xt::xtensor<double, 2>> aos_vals_grady_cart;
        std::vector<xt::xtensor<double, 2>> aos_vals_gradz_cart;
        std::vector<xt::xtensor<double, 2>> aos_vals_gradx_pure;
        std::vector<xt::xtensor<double, 2>> aos_vals_grady_pure;
        std::vector<xt::xtensor<double, 2>> aos_vals_gradz_pure;


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
        static void cDIIS(xt::xtensor<double, 2>& fock, const std::vector<xt::xtensor<double, 2>>& focks, const std::vector<xt::xtensor<double, 2>>& diis_res);
        xt::xtensor<double, 1> compute_rho(const xt::xtensor<double, 2>& Puv, const xt::xtensor<double, 2>& aos_vals);
        std::array<xt::xtensor<double, 1>, 3> compute_rho_grad(const xt::xtensor<double, 2>& Puv, 
                                                               const xt::xtensor<double, 2>& aos_vals, 
                                                               const xt::xtensor<double, 2>& aos_vals_gradx, 
                                                               const xt::xtensor<double, 2>& aos_vals_grady, 
                                                               const xt::xtensor<double, 2>& aos_vals_gradz);
        xt::xtensor<double, 2> xc_quadrature_lda(const xt::xtensor<double, 1>& mweighted_prop, const xt::xtensor<double, 2>& aos_vals);
        xt::xtensor<double, 2> xc_quadrature_gga(const xt::xtensor<double, 1>& mweighted_vsigma, const xt::xtensor<double, 2>& aos_vals, 
                                                 const xt::xtensor<double, 2>& aos_vals_gradx, 
                                                 const xt::xtensor<double, 2>& aos_vals_grady, 
                                                 const xt::xtensor<double, 2>& aos_vals_gradz, 
                                                 const xt::xtensor<double, 1>& rho_gradx, 
                                                 const xt::xtensor<double, 1>& rho_grady, 
                                                 const xt::xtensor<double, 1>& rho_gradz);
        double energy_decomposition(const xt::xtensor<double, 2>& Puv, const xt::xtensor<double, 2>& op);
    }; // end struct DFT


} // end namespace

