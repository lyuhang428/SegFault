#include "../include/hf.hpp"


void sf::HF::cDIIS(xt::xtensor<double, 2>& fock, const std::vector<xt::xtensor<double, 2>>& focks, const std::vector<xt::xtensor<double, 2>>& diis_res)
{
    std::cout << "cDIIS enabled\n";
    assert(fock.shape(0) == fock.shape(1));
    const size_t n = focks.size() + 1;
    const size_t nbf = fock.shape(0);
    xt::xtensor<double, 2> B = xt::zeros<double>({n, n});
    xt::row(B, n-1) = -1.;
    xt::col(B, n-1) = -1.;
    B(n-1, n-1) = 0.;

    for (auto ii=0; ii < n-1; ++ii) {
        for (auto jj=ii; jj < n-1; ++jj) {
            xt::xtensor<double, 2> tmp = xt::zeros<double>({nbf, nbf});
            cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                        CblasTrans, nbf, nbf, 
                        nbf, 1., diis_res[ii].data(), 
                        nbf, diis_res[jj].data(), nbf, 
                        0., tmp.data(), nbf); // diis_res[ii] @ diis_res[jj].T = tmp
            if (jj == ii) {B(ii, jj) = ltr(tmp.data(), nbf); continue;}
            B(ii, jj) = ltr(tmp.data(), nbf);
            B(jj, ii) = B(ii, jj);
        }
    }
    
    xt::xtensor<double, 1> diis_rhs = xt::zeros<double>({n});
    diis_rhs[n-1] = -1.;

    std::vector<int> ipiv(n, 0);
    int info = LAPACKE_dsysv(LAPACK_ROW_MAJOR, 'U', n, 
                            1, B.data(), n, 
                            ipiv.data(), diis_rhs.data(), 1); // B x diis_coef = diis_rhs ; B destroyed, diis_rhs becomes diis_coef
    assert(info == 0);
    fock *= 0.; // set to 0
    for (auto iii=0; iii < n-1; ++iii) fock += focks[iii] * diis_rhs[iii];
    fock = (fock + xt::transpose(fock)) * 0.5;    
}


sf::HF::HF::HF(const std::string& xyzfile, const std::string& name) : xyzfile(xyzfile), name(name)
{
    this->mol = sf::Molecule{xyzfile, name};
}


void sf::HF::HF::scf(int maxiter, double e_convergence, double d_convergence, int nbuffer, const std::string& initial_guess, const bool pure)
{
    const std::vector<libint2::Atom>& atoms = this->mol.atoms;
    const std::vector<libint2::Shell>& shells = pure ? this->mol.shells_pure : this->mol.shells_cart;
    const size_t nbf = pure ? this->mol.nbf_pure : this->mol.nbf_cart;
    
    const int nshell = this->mol.nshell;
    const double e_nuc = this->mol.e_nuc;

    std::vector<size_t> shell2bf = sf::get_shell2bf(shells);
    const size_t       max_nprim = sf::get_max_nprim(shells);
    const int               lmax = sf::get_lmax(shells);

    const xt::xtensor<double, 2> sij = sf::get_olp(shells);
    const xt::xtensor<double, 2> tij = sf::get_kin(shells);
    const xt::xtensor<double, 2> vij = sf::get_ext(shells, atoms);
    const xt::xtensor<double, 2> hij = tij + vij;


    // 重叠矩阵正交化
    // Sij @ u = u @ s
    xt::xtensor<double, 2> u = xt::zeros<double>({nbf, nbf}); // vecs of Sij
    xt::xtensor<double, 1> s = xt::zeros<double>({nbf});      // vals of Sij
    {
        xt::xtensor<double, 2> sij_mutable{sij};                  // sij is const, make a deep copy    
        double sfmin = LAPACKE_dlamch('S'); // safe minium for DSYEVR
        int m;
        std::vector<int> isuppz(2*nbf, 0);
        int info = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, 'V', 'I', 'U', 
                                nbf, sij_mutable.data(), nbf, 0., 
                                0., 1, nbf, 
                                sfmin, &m, s.data(), u.data(), 
                                nbf, isuppz.data()); // sij -> vecs: u ; vals: s
        assert(info == 0);
    }

    xt::xtensor<double, 2> inv_s = xt::zeros<double>({nbf, nbf}); // diag(s^-1/2)
    for (auto i=0; i < nbf; ++i) inv_s(i,i) = std::sqrt(1. / s[i]);
    xt::xtensor<double, 2> sij_inv_half = xt::zeros<double>({nbf, nbf}); // sij_inv_half = (u @ inv_s) @ u.T
    {
        xt::xtensor<double, 2> tmp = xt::zeros<double>({nbf, nbf}); // tmp is local
        cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                    CblasNoTrans, nbf, nbf, 
                    nbf, 1., u.data(), 
                    nbf, inv_s.data(), nbf, 
                    0., tmp.data(), nbf);
        
        cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                    CblasTrans, nbf, nbf, 
                    nbf, 1., tmp.data(), 
                    nbf, u.data(), nbf, 
                    0., sij_inv_half.data(), nbf);
    }

    std::vector<xt::xtensor<double, 2>> focks   ; focks.reserve(maxiter+1);
    std::vector<xt::xtensor<double, 2>> diis_res; diis_res.reserve(maxiter+1);

    int counter = 0;
    double etot_old = std::nan("1");
    xt::xtensor<double, 2> Puv_old = xt::ones<double>({nbf, nbf}) * std::nan("1");
    const int ne = this->mol.ne;
    const int nocc = this->mol.nocc;

    // initial guess Hcore
    xt::xtensor<double, 2> fock = hij;


    /////////////
    //>! SCF <!//
    /////////////
    std::cout << "bf type " << (pure ? "pure" : "cartesian") << std::endl;
    std::cout << std::format("Fock shape = ({:>d}, {:>d})\n", fock.shape(0), fock.shape(1));
    
    while (true) {
        std::cout << "\n>!STEP " << counter+1 << std::endl;

        xt::xtensor<double, 2> fprime = xt::zeros<double>({nbf, nbf}); // fprime = (sij_inv_half @ fock) @ sij_inv_half
        xt::xtensor<double, 2> cprime = xt::zeros<double>({nbf, nbf});
        xt::xtensor<double, 1> e      = xt::zeros<double>({nbf});
        xt::xtensor<double, 2> vecs   = xt::zeros<double>({nbf, nbf});
        
        {
            xt::xtensor<double, 2> tmp = xt::zeros<double>({nbf, nbf});
            cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                        CblasNoTrans, nbf, nbf, 
                        nbf, 1., sij_inv_half.data(), 
                        nbf, fock.data(), nbf, 
                        0., tmp.data(), nbf); // sij_inv_half @ hij = tmp
            cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                        CblasNoTrans, nbf, nbf, 
                        nbf, 1., tmp.data(), 
                        nbf, sij_inv_half.data(), nbf, 
                        0., fprime.data(), nbf); // fprime = (sij_inv_half @ hij) @ sij_inv_half

            double sfmin = LAPACKE_dlamch('S'); // safe minium for DSYEVR
            int m;
            std::vector<int> isuppz(2*nbf, 0);
            int info = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, 'V', 'I', 'U', 
                                    nbf, fprime.data(), nbf, 0., 
                                    0., 1, nbf, 
                                    sfmin, &m, e.data(), cprime.data(), 
                                    nbf, isuppz.data()); // fprime -> vecs: cprime  vals: e
            assert(info == 0);
            cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                        CblasNoTrans, nbf, nbf, 
                        nbf, 1., sij_inv_half.data(), 
                        nbf, cprime.data(), nbf, 
                        0., vecs.data(), nbf); // vecs = sij_inv_half @ cprime
        }

        xt::xtensor<double, 2> Puv = xt::zeros<double>({nbf, nbf}); // Puv = vecs[:,:nocc] @ vecs[:,:nocc].T
        cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                    CblasTrans, nbf, nbf, 
                    nocc, 2., vecs.data(), 
                    nbf, vecs.data(), nbf, 
                    0., Puv.data(), nbf);
        
        auto [Juv, Kuv] = sf::get_JK(shells, Puv); // 已乘密度矩阵

        fock = hij + Juv - Kuv; fock = (fock + xt::transpose(fock)) * 0.5;
        focks.emplace_back(fock);

        // 收敛加速
        // diis_res = sij_inv_half @ ((fock @ Puv) @ sij - sij @ (Puv @ fock)) @ sij_inv_half
        xt::xtensor<double, 2> tmp_lhs = xt::zeros<double>({nbf, nbf});
        xt::xtensor<double, 2> tmp_rhs = xt::zeros<double>({nbf, nbf});
        xt::xtensor<double, 2> lhs     = xt::zeros<double>({nbf, nbf});
        xt::xtensor<double, 2> rhs     = xt::zeros<double>({nbf, nbf});
        xt::xtensor<double, 2> middle  = xt::zeros<double>({nbf, nbf});
        xt::xtensor<double, 2> tmp     = xt::zeros<double>({nbf, nbf});
        xt::xtensor<double, 2> residue = xt::zeros<double>({nbf, nbf});

        cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                    CblasNoTrans, nbf, nbf, 
                    nbf, 1., fock.data(), 
                    nbf, Puv.data(), nbf, 
                    0., tmp_lhs.data(), nbf); // fock @ Puv = tmp_lhs
        cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                    CblasNoTrans, nbf, nbf, 
                    nbf, 1., tmp_lhs.data(), 
                    nbf, sij.data(), nbf, 
                    0., lhs.data(), nbf); // tmp_lhs @ sij = lhs
        cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                    CblasNoTrans, nbf, nbf, 
                    nbf, 1., Puv.data(), 
                    nbf, fock.data(), nbf, 
                    0., tmp_rhs.data(), nbf); // Puv @ fock = tmp_rhs
        cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                    CblasNoTrans, nbf, nbf, 
                    nbf, 1., sij.data(), 
                    nbf, tmp_rhs.data(), nbf, 
                    0., rhs.data(), nbf); // tmp_rhs @ sij = rhs
        middle = lhs - rhs;
        cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                    CblasNoTrans, nbf, nbf, 
                    nbf, 1., sij_inv_half.data(), 
                    nbf, middle.data(), nbf, 
                    0., tmp.data(), nbf); // sij_inv_half @ middle = tmp
        cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                    CblasNoTrans, nbf, nbf, 
                    nbf, 1., tmp.data(), 
                    nbf, sij_inv_half.data(), nbf, 
                    0., residue.data(), nbf); // tmp @ sij_inv_half = residue
    
        diis_res.emplace_back(residue);

        if (counter > 2) cDIIS(fock, focks, diis_res);

        
        // energy decomposition
        xt::xtensor<double, 2> Puv_times_tijT = xt::zeros<double>({nbf, nbf});
        xt::xtensor<double, 2> Puv_times_vijT = xt::zeros<double>({nbf, nbf});
        xt::xtensor<double, 2> Puv_times_JuvT = xt::zeros<double>({nbf, nbf});
        xt::xtensor<double, 2> Puv_times_KuvT = xt::zeros<double>({nbf, nbf});

        cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                    CblasTrans, nbf, nbf, 
                    nbf, 1., Puv.data(), 
                    nbf, tij.data(),  nbf, 
                    0., Puv_times_tijT.data(),  nbf); // Puv @ tij.T
        cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                    CblasTrans, nbf, nbf, 
                    nbf, 1., Puv.data(), 
                    nbf, vij.data(),  nbf, 
                    0., Puv_times_vijT.data(),  nbf); // Puv @ vij,T
        cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                    CblasTrans, nbf, nbf, 
                    nbf, 1., Puv.data(), 
                    nbf, Juv.data(),  nbf, 
                    0., Puv_times_JuvT.data(),  nbf); // Puv @ Juv.T
        cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                    CblasTrans, nbf, nbf, 
                    nbf, 1., Puv.data(), 
                    nbf, Kuv.data(), nbf, 
                    0., Puv_times_KuvT.data(), nbf); // Puv @ Exuv.T
    
        double kin_e         = ltr(Puv_times_tijT.data(),  nbf);
        double ext_e         = ltr(Puv_times_vijT.data(),  nbf);
        double hartree_e     = ltr(Puv_times_JuvT.data(),  nbf) * 0.5;
        double exchange_e    = ltr(Puv_times_KuvT.data(), nbf) * -0.5;
        double etot = 0.5 * xt::sum(Puv * (hij + fock))() + e_nuc;
        double e_diff = etot - etot_old;
        double d_diff = xt::sum(xt::square(Puv - Puv_old))();
        
        std::cout << std::format("One-Electron energy = {:>25.15f}\n", kin_e + ext_e);
        std::cout << std::format("Two-Electron energy = {:>25.15f}\n", hartree_e + exchange_e);
        std::cout << std::format("Total energy        = {:>25.15f}\n", etot);
        std::cout << std::format("Energy difference   = {:>25.15f}\n", e_diff);
        std::cout << std::format("Density difference  = {:>25.15f}\n", d_diff);


        // DIIS
        if (counter > 2) cDIIS(fock, focks, diis_res);


        // convergence check
        if (std::abs(e_diff) <= e_convergence && std::abs(d_diff) <= d_convergence) {
            std::cout << std::endl;
            std::cout << "\nDoubly occupied:\n";
            int iorb = 1;
            for (auto val : xt::view(e, xt::range(0, nocc))) {
                std::cout << std::format("{:>15.10f}    ", val);
                if (iorb % 4 == 0) std::cout << std::endl;
                iorb += 1;
            }

            iorb = 1;
            std::cout << "\nVirtual:\n";
            for (auto val : xt::view(e, xt::range(nocc, nbf))) {
                std::cout << std::format("{:>15.10f}    ", val);
                if (iorb % 4 == 0) std::cout << std::endl;
                iorb += 1;
            }

            std::cout << "\n>! NORMAL TERMINALTION\n";
            break;
        }

        counter++;
        etot_old = etot;
        Puv_old = Puv;

        if (counter > maxiter) {std::cout << "MAXITER exceeded\n"; break;}
    } // end while


    std::cout << ">! safe here\n";
}


