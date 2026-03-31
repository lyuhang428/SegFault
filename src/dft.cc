#include "../include/dft.hpp"



sf::DFT::DFT::DFT(const std::string& xyzfile, const std::string& name) : xyzfile(xyzfile), name(name)
{
    this->mol = sf::Molecule{xyzfile, name};
}


void sf::DFT::DFT::_init_lda()
{
    const size_t nbf_cart = this->nbf_cart;
    const size_t nbf_pure = this->nbf_pure;
    const size_t natom = this->natom;
    const size_t natgrid = this->natgrid;
    const auto grid_global = this->becke.build_grid2();

    this->aos_vals_cart = std::vector<xt::xtensor<double, 2>>{natom, xt::zeros<double>({nbf_cart, natgrid})};
    this->aos_vals_pure = std::vector<xt::xtensor<double, 2>>{natom, xt::zeros<double>({nbf_pure, natgrid})};

    sf::BFs bfs{this->mol};
    const xt::xtensor<double, 2> TF = this->mol.make_tf();

    for (auto iatom=0; iatom < natom; ++iatom) {
        this->aos_vals_cart[iatom] = bfs.get_ao_val(xt::col(grid_global[iatom], 0), xt::col(grid_global[iatom], 1), xt::col(grid_global[iatom], 2));
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf_pure, natgrid, nbf_cart, 1., TF.data(), nbf_cart, this->aos_vals_cart[iatom].data(), natgrid, 0., this->aos_vals_pure[iatom].data(), natgrid);
    }
}


void sf::DFT::DFT::_init_gga()
{
    const size_t nbf_cart = this->nbf_cart;
    const size_t nbf_pure = this->nbf_pure;
    const size_t natom = this->natom;
    const size_t natgrid = this->natgrid;
    const auto grid_global = this->becke.build_grid2();

    this->aos_vals_cart       = std::vector<xt::xtensor<double, 2>>{natom, xt::zeros<double>({nbf_cart, natgrid})};
    this->aos_vals_pure       = std::vector<xt::xtensor<double, 2>>{natom, xt::zeros<double>({nbf_pure, natgrid})};
    this->aos_vals_gradx_cart = std::vector<xt::xtensor<double, 2>>{natom, xt::zeros<double>({nbf_cart, natgrid})};
    this->aos_vals_grady_cart = std::vector<xt::xtensor<double, 2>>{natom, xt::zeros<double>({nbf_cart, natgrid})};
    this->aos_vals_gradz_cart = std::vector<xt::xtensor<double, 2>>{natom, xt::zeros<double>({nbf_cart, natgrid})};
    this->aos_vals_gradx_pure = std::vector<xt::xtensor<double, 2>>{natom, xt::zeros<double>({nbf_pure, natgrid})};
    this->aos_vals_grady_pure = std::vector<xt::xtensor<double, 2>>{natom, xt::zeros<double>({nbf_pure, natgrid})};
    this->aos_vals_gradz_pure = std::vector<xt::xtensor<double, 2>>{natom, xt::zeros<double>({nbf_pure, natgrid})};

    sf::BFs bfs{this->mol};
    const xt::xtensor<double, 2> TF = this->mol.make_tf();

#pragma omp for
    for (auto iatom=0; iatom < natom; ++iatom) {
        this->aos_vals_cart[iatom] = bfs.get_ao_val(xt::col(grid_global[iatom], 0), xt::col(grid_global[iatom], 1), xt::col(grid_global[iatom], 2));
        auto aos_vals_grads_cart = bfs.get_ao_grad(xt::col(grid_global[iatom], 0), xt::col(grid_global[iatom], 1), xt::col(grid_global[iatom], 2));
        this->aos_vals_gradx_cart[iatom] = std::move(aos_vals_grads_cart[0]);
        this->aos_vals_grady_cart[iatom] = std::move(aos_vals_grads_cart[1]);
        this->aos_vals_gradz_cart[iatom] = std::move(aos_vals_grads_cart[2]);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf_pure, natgrid, nbf_cart, 1., TF.data(), nbf_cart, this->aos_vals_cart[iatom].data(), natgrid, 0., this->aos_vals_pure[iatom].data(), natgrid);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf_pure, natgrid, nbf_cart, 1., TF.data(), nbf_cart, this->aos_vals_gradx_cart[iatom].data(), natgrid, 0., this->aos_vals_gradx_pure[iatom].data(), natgrid);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf_pure, natgrid, nbf_cart, 1., TF.data(), nbf_cart, this->aos_vals_grady_cart[iatom].data(), natgrid, 0., this->aos_vals_grady_pure[iatom].data(), natgrid);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf_pure, natgrid, nbf_cart, 1., TF.data(), nbf_cart, this->aos_vals_gradz_cart[iatom].data(), natgrid, 0., this->aos_vals_gradz_pure[iatom].data(), natgrid);
    }
}


void sf::DFT::DFT::_init(int radial_points, int angular_level, int k, bool biased, std::string radial_scheme, int xc_family)
{
    this->becke    = beckegrid::BeckeFuzzyCell{mol.symbols, mol.rms, mol.numbers, mol.xyz, static_cast<size_t>(radial_points), static_cast<size_t>(angular_level),static_cast<size_t>(k), biased};
    this->natom    = mol.natom;
    this->nbf_cart = mol.nbf_cart;
    this->nbf_pure = mol.nbf_pure;
    this->ne       = mol.ne;
    this->nocc     = mol.nocc;

    auto it = std::ranges::find(LEBEDEV_ORDER, angular_level);
    assert(it != std::end(LEBEDEV_ORDER));
    int index = it - std::begin(LEBEDEV_ORDER);

    this->nang    = LEBEDEV_LEVEL[index];
    this->nrad    = static_cast<int>(becke.nrad);
    this->natgrid = this->nrad * this->nang;
    this->ngrid   = this->natom * this->nrad * this->nang;
    
    if      (xc_family == XC_FAMILY_LDA) {std::cout << "Initialize LDA functional\n"; this->_init_lda();}
    else if (xc_family == XC_FAMILY_GGA) {std::cout << "Initialize GGA functional\n"; this->_init_gga();}
    else {std::cerr << "Unknown or un-supported xc\n"; exit(-1);}
}


xt::xtensor<double, 1> sf::DFT::DFT::compute_rho(const xt::xtensor<double, 2>& Puv, const xt::xtensor<double, 2>& aos_vals)
{
    const size_t nbf = Puv.shape(0);
    const size_t natgrid = aos_vals.shape(1);
    xt::xtensor<double, 2> tmp = xt::zeros<double>({nbf, natgrid});
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, natgrid, nbf, 1., Puv.data(), nbf, aos_vals.data(), natgrid, 0., tmp.data(), natgrid);
    return xt::sum(aos_vals * tmp, {0});
}


std::array<xt::xtensor<double, 1>, 3> sf::DFT::DFT::compute_rho_grad(const xt::xtensor<double, 2>& Puv, 
                                                                     const xt::xtensor<double, 2>& aos_vals, 
                                                                     const xt::xtensor<double, 2>& aos_vals_gradx, 
                                                                     const xt::xtensor<double, 2>& aos_vals_grady, 
                                                                     const xt::xtensor<double, 2>& aos_vals_gradz)
{
    const size_t nbf = Puv.shape(0);
    const size_t natgrid = aos_vals.shape(1);
    xt::xtensor<double, 2> tmp_x = xt::zeros<double>({nbf, natgrid});
    xt::xtensor<double, 2> tmp_y = xt::zeros<double>({nbf, natgrid});
    xt::xtensor<double, 2> tmp_z = xt::zeros<double>({nbf, natgrid});

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, natgrid, nbf, 1., Puv.data(), nbf, aos_vals_gradx.data(), natgrid, 0., tmp_x.data(), natgrid);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, natgrid, nbf, 1., Puv.data(), nbf, aos_vals_grady.data(), natgrid, 0., tmp_y.data(), natgrid);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, natgrid, nbf, 1., Puv.data(), nbf, aos_vals_gradz.data(), natgrid, 0., tmp_z.data(), natgrid);
    
    xt::xtensor<double, 1> rho_gradx = 2. * xt::sum(aos_vals * tmp_x, {0});
    xt::xtensor<double, 1> rho_grady = 2. * xt::sum(aos_vals * tmp_y, {0});
    xt::xtensor<double, 1> rho_gradz = 2. * xt::sum(aos_vals * tmp_z, {0});

    return {rho_gradx, rho_grady, rho_gradz};
}



xt::xtensor<double, 2> sf::DFT::DFT::xc_quadrature_lda(const xt::xtensor<double, 1>& mweighted_prop, const xt::xtensor<double, 2>& aos_vals)
{
    const size_t nbf = aos_vals.shape(0);
    const size_t natgrid = aos_vals.shape(1);
    xt::xtensor<double, 2> mweighted_prop_aos_vals = aos_vals * mweighted_prop;
    xt::xtensor<double, 2> mat_local = xt::zeros<double>({nbf, nbf});
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nbf, nbf, natgrid, 1., mweighted_prop_aos_vals.data(), natgrid, aos_vals.data(), natgrid, 0., mat_local.data(), nbf);
    return mat_local;
}


xt::xtensor<double, 2> sf::DFT::DFT::xc_quadrature_gga(const xt::xtensor<double, 1>& mweighted_vsigma, 
                                                       const xt::xtensor<double, 2>& aos_vals, 
                                                       const xt::xtensor<double, 2>& aos_vals_gradx, 
                                                       const xt::xtensor<double, 2>& aos_vals_grady, 
                                                       const xt::xtensor<double, 2>& aos_vals_gradz, 
                                                       const xt::xtensor<double, 1>& rho_gradx, 
                                                       const xt::xtensor<double, 1>& rho_grady, 
                                                       const xt::xtensor<double, 1>& rho_gradz)
{
    const size_t nbf = aos_vals.shape(0);
    const size_t natgrid = aos_vals.shape(1);
    xt::xtensor<double, 2> mweighted_vsigma_aos_vals = aos_vals * mweighted_vsigma;
    xt::xtensor<double, 2> aos_vals_grad_dot_rho_grad = aos_vals_gradx * rho_gradx + aos_vals_grady * rho_grady + aos_vals_gradz * rho_gradz;
    xt::xtensor<double, 2> mat_local = xt::zeros<double>({nbf, nbf});
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nbf, nbf, natgrid, 1., mweighted_vsigma_aos_vals.data(), natgrid, aos_vals_grad_dot_rho_grad.data(), natgrid, 0., mat_local.data(), nbf);
    return 2. * mat_local;
}



double sf::DFT::DFT::energy_decomposition(const xt::xtensor<double, 2>& Puv, const xt::xtensor<double, 2>& op)
{
    const size_t nbf = Puv.shape(0);
    xt::xtensor<double, 2> tmp = xt::zeros<double>({nbf, nbf});
    // Puv @ Op.T
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nbf, nbf, nbf, 1., Puv.data(), nbf, op.data(), nbf, 0., tmp.data(), nbf);
    return ltr(tmp.data(), nbf);
}


void sf::DFT::DFT::cDIIS(xt::xtensor<double, 2>& fock, const std::vector<xt::xtensor<double, 2>>& focks, const std::vector<xt::xtensor<double, 2>>& diis_res)
{
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
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nbf, nbf, nbf, 1., diis_res[ii].data(), nbf, diis_res[jj].data(), nbf, 0., tmp.data(), nbf);
            if (jj == ii) {B(ii, jj) = ltr(tmp.data(), nbf); continue;}
            B(ii, jj) = ltr(tmp.data(), nbf);
            B(jj, ii) = B(ii, jj);
        }
    }
    
    xt::xtensor<double, 1> diis_rhs = xt::zeros<double>({n});
    diis_rhs[n-1] = -1.;

    std::vector<int> ipiv(n, 0);
    int info = LAPACKE_dsysv(LAPACK_ROW_MAJOR, 'U', n, 1, B.data(), n, ipiv.data(), diis_rhs.data(), 1);
    assert(info == 0);
    fock *= 0.;
    for (auto iii=0; iii < n-1; ++iii) fock += focks[iii] * diis_rhs[iii];
    fock = (fock + xt::transpose(fock)) * 0.5;    
}


void sf::DFT::DFT::scf(int maxiter, double e_convergence, double d_convergence, 
                       int nbuffer, std::string initial_guess, 
                       int X_id, int C_id, bool pure, 
                       int radial_points, int angular_level, 
                       int k, bool biased, std::string radial_scheme)
{
    int family_x, family_c, number_x, number_c;
    int kind_x = xc_family_from_id(X_id, &family_x, &number_x);
    int kind_c = xc_family_from_id(C_id, &family_c, &number_c);
    assert(family_x == family_c);
    
#ifdef TIMEIT
const auto t0_init{std::chrono::steady_clock::now()};
#endif
    this->_init(radial_points, angular_level, k, biased, "becke", family_x);
#ifdef TIMEIT
const auto t1_init{std::chrono::steady_clock::now()};
std::cout << "init      " << std::chrono::duration<double>(t1_init - t0_init) << std::endl;
#endif

    const size_t nbf = pure ? this->nbf_pure : this->nbf_cart;
    const size_t natom = this->natom;
    const size_t natgrid = this->natgrid;
    const size_t ngrid = this->ngrid;
    const std::vector<libint2::Atom>& atoms = this->mol.atoms;
    const std::vector<xt::xtensor<double, 2>>& aos_vals       = pure ? this->aos_vals_pure : this->aos_vals_cart;
    const std::vector<xt::xtensor<double, 2>>& aos_vals_gradx = pure ? this->aos_vals_gradx_pure : this->aos_vals_gradx_cart;
    const std::vector<xt::xtensor<double, 2>>& aos_vals_grady = pure ? this->aos_vals_grady_pure : this->aos_vals_grady_cart;
    const std::vector<xt::xtensor<double, 2>>& aos_vals_gradz = pure ? this->aos_vals_gradz_pure : this->aos_vals_gradz_cart;

    const xt::xtensor<double, 2> sij = pure ? sf::get_olp(this->mol.shells_pure)        : sf::get_olp(this->mol.shells_cart);
    const xt::xtensor<double, 2> tij = pure ? sf::get_kin(this->mol.shells_pure)        : sf::get_kin(this->mol.shells_cart);
    const xt::xtensor<double, 2> vij = pure ? sf::get_ext(this->mol.shells_pure, atoms) : sf::get_ext(this->mol.shells_cart, atoms);
    const xt::xtensor<double, 2> hij = tij + vij;
    assert(sij.shape(0) == nbf && sij.shape(1) == nbf);

    xt::xtensor<double, 2> mweights = xt::zeros<double>({natom, natgrid});
    for (auto iatom=0; iatom < natom; ++iatom) {
        xt::xtensor<double, 1> wrad = xt::square(xt::row(this->becke.xrwcheb[iatom], 1)) * xt::row(this->becke.xrwcheb[iatom], 2);
        xt::xtensor<double, 1> wang = xt::row(this->becke.xwleb, 3);
        xt::xtensor<double, 2> outer = xt::zeros<double>({wrad.size(), wang.size()});
        cblas_dger(CblasRowMajor, wrad.size(), wang.size(), 1., wrad.data(), 1, wang.data(), 1, outer.data(), outer.shape(1));
        xt::row(mweights, iatom) = xt::ravel(outer) * this->becke.weights[iatom];
    }

    this->header_log();

    int counter = 0;
    double etot_old = std::nan("1");
    xt::xtensor<double, 2> Puv_old = xt::ones<double>({nbf, nbf}) * std::nan("1");
    std::vector<xt::xtensor<double, 2>>    focks;    focks.reserve(maxiter+1);
    std::vector<xt::xtensor<double, 2>> diis_res; diis_res.reserve(maxiter+1);


    xt::xtensor<double, 2> u = xt::zeros<double>({nbf, nbf});
    xt::xtensor<double, 1> s = xt::zeros<double>({nbf});
    
    {
        xt::xtensor<double, 2> sij_mutable{sij};
        double sfmin = LAPACKE_dlamch('S');
        int m;
        std::vector<int> isuppz(2*nbf, 0);
        int info = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, 'V', 'I', 'U', nbf, sij_mutable.data(), nbf, 0., 0., 1, nbf, sfmin, &m, s.data(), u.data(), nbf, isuppz.data());
        assert(info == 0);
    }

    xt::xtensor<double, 2> inv_s = xt::zeros<double>({nbf, nbf});
    for (auto i=0; i < nbf; ++i) inv_s(i,i) = std::sqrt(1. / s[i]);
    xt::xtensor<double, 2> sij_inv_half = xt::zeros<double>({nbf, nbf}); 
    
    {
        xt::xtensor<double, 2> tmp = xt::zeros<double>({nbf, nbf});
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, nbf, nbf, 1., u.data(), nbf, inv_s.data(), nbf, 0., tmp.data(), nbf);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nbf, nbf, nbf, 1., tmp.data(), nbf, u.data(), nbf, 0., sij_inv_half.data(), nbf);
    }

 
    xt::xtensor<double, 2> fprime = xt::zeros<double>({nbf, nbf});
    xt::xtensor<double, 2> cprime = xt::zeros<double>({nbf, nbf});
    xt::xtensor<double, 1> e = xt::zeros<double>({nbf});
    xt::xtensor<double, 2> vecs = xt::zeros<double>({nbf, nbf});
    if (initial_guess == "core") {
        xt::xtensor<double, 2> tmp = xt::zeros<double>({nbf, nbf});
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, nbf, nbf, 1., sij_inv_half.data(), nbf, hij.data(), nbf, 0., tmp.data(), nbf);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, nbf, nbf, 1., tmp.data(), nbf, sij_inv_half.data(), nbf, 0., fprime.data(), nbf);

        xt::xtensor<double, 2> fprime_mutable{fprime};
        double sfmin = LAPACKE_dlamch('S');
        int m;
        std::vector<int> isuppz(2*nbf, 0);
        int info = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, 'V', 'I', 'U', nbf, fprime_mutable.data(), nbf, 0., 0., 1, nbf, sfmin, &m, e.data(), cprime.data(), nbf, isuppz.data());
        assert(info == 0);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, nbf, nbf, 1., sij_inv_half.data(), nbf, cprime.data(), nbf, 0., vecs.data(), nbf);
    }
    else if(initial_guess == "SAD") {std::cerr << "SAD not implemented yet\n"; exit(-1);}
    else {std::cerr << "Unknown initial guess\n"; exit(-1);}



    xt::xtensor<double, 2> Puv = xt::zeros<double>({nbf, nbf});
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nbf, nbf, this->nocc, 2., vecs.data(), nbf, vecs.data(), nbf, 0., Puv.data(), nbf);


    xc_func_type funcx, funcc;
    xc_func_init(&funcx, X_id, XC_UNPOLARIZED);
    xc_func_init(&funcc, C_id, XC_UNPOLARIZED);
    

    
    std::cout << std::string(65, '=') << std::endl;
    std::cout << std::string(30, '=') << " SCF " << std::string(30, '=') << std::endl;
    std::cout << std::string(65, '=') << std::endl;
    std::cout << "bf type " << (pure ? "pure" : "cartesian") << std::endl;
    std::cout << "Fock shape " << "(" << hij.shape(0) << ", " << hij.shape(1) << ")" << std::endl;
    std::cout << std::endl;
    std::cout << "Step |    1e energy    |    2e energy    |    xc energy    |    total    |    ne    |    ΔE    |    ΔD\n";


#ifdef TIMEIT
    auto t0_eden_grad{std::chrono::steady_clock::now()};
    auto t1_eden_grad{std::chrono::steady_clock::now()};
    auto t0_xc{std::chrono::steady_clock::now()};
    auto t1_xc{std::chrono::steady_clock::now()};
#endif


    while (true) {
    xt::xtensor<double, 2> rho = xt::zeros<double>({natom, natgrid});
    for (auto iatom=0; iatom < natom; ++iatom) xt::row(rho, iatom) = this->compute_rho(Puv, aos_vals[iatom]);
    double ne_quad = xt::sum(rho * mweights)();
    rho *= this->ne / ne_quad;

#ifdef TIMEIT
const auto t0_eri{std::chrono::steady_clock::now()};
#endif
    xt::xtensor<double, 2> Juv = pure ? sf::get_J_par(this->mol.shells_pure, Puv) : sf::get_J_par(this->mol.shells_cart, Puv);    
#ifdef TIMEIT
const auto t1_eri{std::chrono::steady_clock::now()};
#endif
    xt::xtensor<double, 2> Kuv  = xt::zeros<double>({nbf, nbf});
    xt::xtensor<double, 2> Cuv  = xt::zeros<double>({nbf, nbf});
    xt::xtensor<double, 2> Exuv = xt::zeros<double>({nbf, nbf});
    xt::xtensor<double, 2> Ecuv = xt::zeros<double>({nbf, nbf});
    
    xt::xtensor<double, 2> vx = xt::zeros<double>({natom, natgrid});
    xt::xtensor<double, 2> ex = xt::zeros<double>({natom, natgrid});
    xt::xtensor<double, 2> vc = xt::zeros<double>({natom, natgrid});
    xt::xtensor<double, 2> ec = xt::zeros<double>({natom, natgrid});
    
    
    if (family_x == 2) { // GGA
        std::vector<xt::xtensor<double, 1>> rho_gradx{natom, xt::zeros<double>({natgrid})};
        std::vector<xt::xtensor<double, 1>> rho_grady{natom, xt::zeros<double>({natgrid})};
        std::vector<xt::xtensor<double, 1>> rho_gradz{natom, xt::zeros<double>({natgrid})};
        xt::xtensor<double, 2> sigma   = xt::zeros<double>({natom, natgrid});
        xt::xtensor<double, 2> vsigmax = xt::zeros<double>({natom, natgrid});
        xt::xtensor<double, 2> vsigmac = xt::zeros<double>({natom, natgrid});
        
#ifdef TIMEIT
t0_eden_grad = std::chrono::steady_clock::now();
#endif
#pragma omp for
        for (auto iatom=0; iatom < natom; ++iatom) {
            auto rho_grads = this->compute_rho_grad(Puv, aos_vals[iatom], aos_vals_gradx[iatom], aos_vals_grady[iatom], aos_vals_gradz[iatom]);
            rho_gradx[iatom] = rho_grads[0];
            rho_grady[iatom] = rho_grads[1];
            rho_gradz[iatom] = rho_grads[2];
            xt::row(sigma, iatom) = xt::square(rho_grads[0]) + xt::square(rho_grads[1]) + xt::square(rho_grads[2]);
        }
#ifdef TIMEIT
t1_eden_grad = std::chrono::steady_clock::now();
#endif
        
        xc_gga_vxc(&funcx, ngrid, rho.data(), sigma.data(), vx.data(), vsigmax.data());
        xc_gga_exc(&funcx, ngrid, rho.data(), sigma.data(), ex.data());
        xc_gga_vxc(&funcc, ngrid, rho.data(), sigma.data(), vc.data(), vsigmac.data());
        xc_gga_exc(&funcc, ngrid, rho.data(), sigma.data(), ec.data());
        

#ifdef TIMEIT
t0_xc = std::chrono::steady_clock::now();
#endif
#pragma omp parallel for
        for (auto iatom=0; iatom < natom; ++iatom) {
            xt::xtensor<double, 1> vx_mweights      = xt::row(mweights, iatom) * xt::row(vx,      iatom);
            xt::xtensor<double, 1> vc_mweights      = xt::row(mweights, iatom) * xt::row(vc,      iatom);
            xt::xtensor<double, 1> ex_mweights      = xt::row(mweights, iatom) * xt::row(ex,      iatom);
            xt::xtensor<double, 1> ec_mweights      = xt::row(mweights, iatom) * xt::row(ec,      iatom);
            xt::xtensor<double, 1> vsigmax_mweights = xt::row(mweights, iatom) * xt::row(vsigmax, iatom);
            xt::xtensor<double, 1> vsigmac_mweights = xt::row(mweights, iatom) * xt::row(vsigmac, iatom);
            
            xt::xtensor<double, 2> Kuv_local = this->xc_quadrature_lda(vx_mweights, aos_vals[iatom]);
            Kuv_local += 2. * this->xc_quadrature_gga(vsigmax_mweights, aos_vals[iatom], 
                                                      aos_vals_gradx[iatom], aos_vals_grady[iatom], aos_vals_gradz[iatom], 
                                                      rho_gradx[iatom], rho_grady[iatom], rho_gradz[iatom]);
            
            xt::xtensor<double, 2> Cuv_local = this->xc_quadrature_lda(vc_mweights, aos_vals[iatom]);
            Cuv_local += 2. * this->xc_quadrature_gga(vsigmac_mweights, aos_vals[iatom], 
                                                      aos_vals_gradx[iatom], aos_vals_grady[iatom], aos_vals_gradz[iatom], 
                                                      rho_gradx[iatom], rho_grady[iatom], rho_gradz[iatom]);
            
            xt::xtensor<double, 2> Exuv_local = this->xc_quadrature_lda(ex_mweights, aos_vals[iatom]);
            xt::xtensor<double, 2> Ecuv_local = this->xc_quadrature_lda(ec_mweights, aos_vals[iatom]);
#pragma omp critical
            {
            Kuv += Kuv_local;
            Cuv += Cuv_local;
            Exuv += Exuv_local;
            Ecuv += Ecuv_local;
            }
        }
#ifdef TIMEIT
t1_xc = std::chrono::steady_clock::now();
#endif
    }
    else if (family_x == 1) { // LDA
        xc_lda_vxc(&funcx, ngrid, rho.data(), vx.data());
        xc_lda_exc(&funcx, ngrid, rho.data(), ex.data());
        xc_lda_vxc(&funcc, ngrid, rho.data(), vc.data());
        xc_lda_exc(&funcc, ngrid, rho.data(), ec.data());
        
#ifdef TIMEIT
t0_xc = std::chrono::steady_clock::now();
#endif
#pragma omp for
        for (auto iatom=0; iatom < natom; ++iatom) {
            xt::xtensor<double, 1> vx_mweights = xt::row(mweights, iatom) * xt::row(vx, iatom);
            xt::xtensor<double, 1> vc_mweights = xt::row(mweights, iatom) * xt::row(vc, iatom);
            xt::xtensor<double, 1> ex_mweights = xt::row(mweights, iatom) * xt::row(ex, iatom);
            xt::xtensor<double, 1> ec_mweights = xt::row(mweights, iatom) * xt::row(ec, iatom);
            
            xt::xtensor<double, 2> Kuv_local  = this->xc_quadrature_lda(vx_mweights, aos_vals[iatom]);
            xt::xtensor<double, 2> Cuv_local  = this->xc_quadrature_lda(vc_mweights, aos_vals[iatom]);
            xt::xtensor<double, 2> Exuv_local = this->xc_quadrature_lda(ex_mweights, aos_vals[iatom]);
            xt::xtensor<double, 2> Ecuv_local = this->xc_quadrature_lda(ec_mweights, aos_vals[iatom]);
#pragma omp critical
            {
            Kuv  += Kuv_local;
            Cuv  += Cuv_local;
            Exuv += Exuv_local;
            Ecuv += Ecuv_local;
            }
        }
#ifdef TIMEIT
t1_xc = std::chrono::steady_clock::now();
#endif
    }
    else {std::cerr << "Unknown or un-supported xc\n"; exit(-1);}


    xt::xtensor<double, 2> fock = tij + vij + Juv + Kuv + Cuv;
    fock = (fock + xt::transpose(fock)) * 0.5;
    focks.emplace_back(fock);

    double kin_e         = energy_decomposition(Puv, tij);
    double ext_e         = energy_decomposition(Puv, vij);
    double hartree_e     = energy_decomposition(Puv, Juv) * 0.5;
    double exchange_e    = energy_decomposition(Puv, Exuv);
    double correlation_e = energy_decomposition(Puv, Ecuv);
    double etot          = kin_e + ext_e + hartree_e + exchange_e + correlation_e + this->mol.e_nuc;

    std::cout << std::format("{:<4d}   {:>.12f}   {:>.12f}   {:>.12f}   {:>.12f}   {:>.12f}  ", 
    counter+1, kin_e+ext_e, hartree_e, exchange_e+correlation_e, etot, ne_quad);

    xt::xtensor<double, 2> tmp_lhs = xt::zeros<double>({nbf, nbf});
    xt::xtensor<double, 2> tmp_rhs = xt::zeros<double>({nbf, nbf});
    xt::xtensor<double, 2> lhs     = xt::zeros<double>({nbf, nbf});
    xt::xtensor<double, 2> rhs     = xt::zeros<double>({nbf, nbf});
    xt::xtensor<double, 2> middle  = xt::zeros<double>({nbf, nbf});
    xt::xtensor<double, 2> tmp     = xt::zeros<double>({nbf, nbf});
    xt::xtensor<double, 2> residue = xt::zeros<double>({nbf, nbf});

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, nbf, nbf, 1., fock.data(), nbf, Puv.data(), nbf, 0., tmp_lhs.data(), nbf);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, nbf, nbf, 1., tmp_lhs.data(), nbf, sij.data(), nbf, 0., lhs.data(), nbf);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, nbf, nbf, 1., Puv.data(), nbf, fock.data(), nbf, 0., tmp_rhs.data(), nbf);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, nbf, nbf, 1., sij.data(), nbf, tmp_rhs.data(), nbf, 0., rhs.data(), nbf);
    middle = lhs - rhs;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, nbf, nbf, 1., sij_inv_half.data(), nbf, middle.data(), nbf, 0., tmp.data(), nbf); 
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, nbf, nbf, 1., tmp.data(), nbf, sij_inv_half.data(), nbf, 0., residue.data(), nbf); 
    diis_res.emplace_back(residue);
    if (counter >= 2) this->cDIIS(fock, focks, diis_res);


    {
        xt::xtensor<double, 2> tmp = xt::zeros<double>({nbf, nbf});
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, nbf, nbf, 1., sij_inv_half.data(), nbf, fock.data(), nbf, 0., tmp.data(), nbf);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, nbf, nbf, 1., tmp.data(), nbf, sij_inv_half.data(), nbf, 0., fprime.data(), nbf);
        xt::xtensor<double, 2> fprime_mutable{fprime};
        double sfmin = LAPACKE_dlamch('S');
        int m;
        std::vector<int> isuppz(2*nbf, 0);
        int info = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, 'V', 'I', 'U', 
                                nbf, fprime_mutable.data(), nbf, 0., 
                                0., 1, nbf, 
                                sfmin, &m, e.data(), cprime.data(), 
                                nbf, isuppz.data());
        assert(info == 0);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nbf, nbf, nbf, 1., sij_inv_half.data(), nbf, cprime.data(), nbf, 0., vecs.data(), nbf);
    }

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nbf, nbf, this->nocc, 2., vecs.data(), nbf, vecs.data(), nbf, 0., Puv.data(), nbf);

    double e_diff = etot - etot_old;
    double d_diff = xt::sum(xt::square(Puv - Puv_old))();
    std::cout << std::format(" {:>.12f}   {:>.12f} \n", e_diff, d_diff);
    if (e_diff < e_convergence && d_diff < d_convergence) {
        std::cout << std::format("\n!SCF converged after {} step. Total energy = {:<.15f} a.u.\n", counter+1, etot);
        
        int iorb = 1;
        std::cout << "\nDoubly occupied:\n";
        for (auto val : xt::view(e, xt::range(0, this->nocc))) {
            std::cout << std::setw(15) << std::setprecision(10) << std::fixed << std::right << val << "    ";
            if (iorb % 4 == 0) std::cout << std::endl;
            iorb += 1;
        }
        
        iorb = 1;
        std::cout << "\nVirtual:\n";
        for (auto val : xt::view(e, xt::range(this->nocc, nbf))) {
            std::cout << std::setw(15) << std::setprecision(10) << std::fixed << std::right << val << "    ";
            if (iorb % 4 == 0) std::cout << std::endl;
            iorb += 1;
        }
        std::cout << "\n>! JOB DONE" << std::endl;
    
        break;
    }

    if (counter >= maxiter) {
        std::cout << "SCF did not converge in " << maxiter << " step. Energy in the last iteration " << etot << std::endl;
        break;
    }


    counter++;
    Puv_old = Puv;
    etot_old = etot;


#ifdef TIMEIT
std::cout << "eri       " << std::chrono::duration<double>(t1_eri - t0_eri) << std::endl;
if (family_x == 2) std::cout << "eden grad " << std::chrono::duration<double>(t1_eden_grad - t0_eden_grad) << std::endl;
std::cout << "xc        " << std::chrono::duration<double>(t1_xc - t0_xc) << std::endl;
#endif
    } // end while


    xc_func_end(&funcx);
    xc_func_end(&funcc);

    std::cout << "safe here\n";
} // end DFT.scf()

