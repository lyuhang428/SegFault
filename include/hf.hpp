#ifndef INCLUDE_HF_HPP
#define INCLUDE_HF_HPP

#include <iostream>
#include <format>
#include <iomanip>
#include <vector>

#include "xtensor/containers/xtensor.hpp"
#include "xtensor/views/xview.hpp"
#include "xtensor/generators/xgenerator.hpp"
#include "xtensor/io/xnpy.hpp"
#include "libint2.hpp"

#include "mints.hpp"
#include "constants.hpp"
#include "lblas.hpp"

namespace sf::HF {
void cDIIS(xt::xtensor<double, 2>& fock, const std::vector<xt::xtensor<double, 2>>& focks, const std::vector<xt::xtensor<double, 2>>& diis_res);

struct HF {
    HF(const std::string& xyzfile, const std::string& name); //>! mol.xyz base.g94
    HF() = default;
    ~HF() = default;

    void scf(int               maxiter=30,
            double             e_convergence=1e-6,
            double             d_convergence=1e-6,
            int                nbuffer=-1,
            const std::string& initial_guess="core", 
            const bool         pure=true);

private:

    std::string                   xyzfile;
    std::string                      name;
    sf::Molecule                      mol;
    std::vector<double>  kinetic_energies;
    std::vector<double> external_energies;
    std::vector<double>  hartree_energies;
    std::vector<double> exchange_energies;
    std::vector<double>             etots;
};



} // end namespace sf::HF


#endif
