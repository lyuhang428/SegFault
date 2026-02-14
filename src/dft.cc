#include "../include/dft.hpp"

#ifdef TIMEIT
#include <chrono>
#endif



sf::DFT::DFT::DFT(const std::string& xyzfile, const std::string& name) : xyzfile(xyzfile), name(name)
{
    this->mol = sf::Molecule{xyzfile, name};
}


void sf::DFT::DFT::header_log() const
{
    std::cout << std::format("Number of atoms          = {:<d}\n", this->natom);
    std::cout << std::format("Number of electrons      = {:<d}\n", this->ne);
    std::cout << std::format("Number of radial grids   = {:<d}\n", this->nrad);
    std::cout << std::format("Number of angular grids  = {:<d}\n", this->nang);
    std::cout << std::format("Number of total grids    = {:<d}\n", this->ngrid);
    std::cout << std::format("Nuclear repulsion energy = {:<.15f} a.u.\n", this->mol.e_nuc);
    std::cout << std::endl;
}


void sf::DFT::DFT::energy_log() const
{
    std::cout << std::format("Number of electron via quadrature = {:>.12f}\n", this->ne_quads.back());
    std::cout << std::format("One-electron energy               = {:>.12f}\n", this->kinetic_energies.back() + this->external_energies.back());
    std::cout << std::format("Two-electron energy               = {:>.12f}\n", this->hartree_energies.back());
    std::cout << std::format("Exchange-correlation energy       = {:>.12f}    {:>.12f}    {:>.12f}\n", this->exchange_energies.back() + this->correlation_energies.back(), this->exchange_energies.back(), this->correlation_energies.back());
    std::cout << std::format("Total energy                      = {:>.12f}\n", this->etots.back());
    std::cout << std::format("Energy difference                 = {:>.12f}\n", this->etots.back() - this->etots[this->etots.size()-2]);
}


void sf::DFT::DFT::init(const int          radial_points,
                        const int          angular_level,
                        const int          k,
                        const bool         biased,
                        const std::string& radial_scheme)
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
    
    // just too lazy to type static_cast<size_t>
    const size_t nbf_cart = this->nbf_cart;
    const size_t nbf_pure = this->nbf_pure;
    const size_t natom = this->natom;
    const size_t natgrid = this->natgrid;
    
    const auto grid_global = this->becke.build_grid2();
    xt::xtensor<double, 2> grid_total = xt::zeros<double>({this->ngrid, 3});
    for (auto iatom=0; iatom < this->natom; ++iatom)
    {
        auto view_tmp = xt::view(grid_total, xt::range(iatom * this->natgrid, (iatom+1) * this->natgrid), xt::all());
        view_tmp = grid_global[iatom];
    }
    
    // [natom, (nbf, natgrid)]
    this->aos_vals_cart = std::vector<xt::xtensor<double, 2>>{natom, xt::zeros<double>({nbf_cart, natgrid})};
    this->aos_vals_pure = std::vector<xt::xtensor<double, 2>>{natom, xt::zeros<double>({nbf_pure, natgrid})};

    sf::BFs bfs{this->mol};
    const xt::xtensor<double, 2> TF = this->mol.make_tf();

#pragma omp parallel for schedule(dynamic)
    for (auto iatom=0; iatom < this->natom; ++iatom) {
        auto x_view = xt::col(grid_global[iatom], 0); // x
        auto y_view = xt::col(grid_global[iatom], 1); // y
        auto z_view = xt::col(grid_global[iatom], 2); // z
        auto ao_val_tmp = bfs.ao_val(x_view, y_view, z_view); // (nbf, natgrid)
        this->aos_vals_cart[iatom] = ao_val_tmp;
        xt::xtensor<double, 2> view_cart_mat{this->aos_vals_cart[iatom]};
        xt::xtensor<double, 2> view_pure_mat = xt::zeros<double>({TF.shape(0), view_cart_mat.shape(1)});
        cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                    CblasNoTrans, TF.shape(0), view_cart_mat.shape(1), 
                    view_cart_mat.shape(0), 1., TF.data(), 
                    TF.shape(1), view_cart_mat.data(), view_cart_mat.shape(1), 
                    0., view_pure_mat.data(), view_pure_mat.shape(1)); // TF @ view_cart_map = view_pure_mat
        this->aos_vals_pure[iatom] = view_pure_mat;
    }
} // end DFT.init()


void sf::DFT::DFT::cDIIS(xt::xtensor<double, 2>& fock, const std::vector<xt::xtensor<double, 2>>& focks, const std::vector<xt::xtensor<double, 2>>& diis_res)
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


void sf::DFT::DFT::scf(const int maxiter, 
                       const double e_convergence, 
                       const double d_convergence, 
                       const int nbuffer, 
                       const std::string& initial_guess, 
                       const int X_id, 
                       const int C_id, 
                       const bool pure)
{
    const std::vector<libint2::Atom>& atoms = this->mol.atoms;
    const std::vector<xt::xtensor<double, 2>>& aos_vals = pure ? this->aos_vals_pure : this->aos_vals_cart;
    const size_t nbf = pure ? this->nbf_pure : this->nbf_cart;

    const xt::xtensor<double, 2> sij = pure ? sf::get_olp(this->mol.shells_pure)        : sf::get_olp(this->mol.shells_cart);
    const xt::xtensor<double, 2> tij = pure ? sf::get_kin(this->mol.shells_pure)        : sf::get_kin(this->mol.shells_cart);
    const xt::xtensor<double, 2> vij = pure ? sf::get_ext(this->mol.shells_pure, atoms) : sf::get_ext(this->mol.shells_cart, atoms);
    const xt::xtensor<double, 2> hij = tij + vij;
    assert(sij.shape(0) == nbf && sij.shape(1) == nbf);

    xt::xtensor<double, 2> mweights = xt::zeros<double>({this->natom, this->natgrid});
    for (auto iatom=0; iatom < this->natom; ++iatom) {
        // (nrad, ) 径向权重
        // (nang, ) 角度权重
        // outer (nrad, nang) 外积
        xt::xtensor<double, 1> wrad = xt::square(xt::row(this->becke.xrwcheb[iatom], 1)) * xt::row(this->becke.xrwcheb[iatom], 2);
        xt::xtensor<double, 1> wang = xt::row(this->becke.xwleb, 3);
        xt::xtensor<double, 2> outer = xt::zeros<double>({wrad.size(), wang.size()});
        cblas_dger(CblasRowMajor, wrad.size(), wang.size(), 
                    1., wrad.data(), 1, 
                    wang.data(), 1, outer.data(), outer.shape(1)); // np.outer(wrad, wang)
        xt::row(mweights, iatom) = xt::ravel(outer) * this->becke.weights[iatom];
    }
    

    this->header_log();

    
    // 参数设置
    int counter = 0;
    this->etots.emplace_back(std::nan("1"));
    xt::xtensor<double, 2> Puv_old = xt::ones<double>({nbf, nbf}) * std::nan("1");
    std::vector<xt::xtensor<double, 2>> focks;
    std::vector<xt::xtensor<double, 2>> diis_res;
    focks.reserve(maxiter+1);
    diis_res.reserve(maxiter+1);
    this->kinetic_energies.reserve(maxiter + 1);
    this->external_energies.reserve(maxiter + 1);
    this->hartree_energies.reserve(maxiter + 1);
    this->exchange_energies.reserve(maxiter + 1);
    this->correlation_energies.reserve(maxiter + 1);
    this->etots.reserve(maxiter + 1);
    this->ne_quads.reserve(maxiter + 1);


    // 重叠矩阵正交化
    // Sij @ u = u @ s
    xt::xtensor<double, 2> u = xt::zeros<double>({nbf, nbf}); // vecs of Sij
    xt::xtensor<double, 1> s = xt::zeros<double>({nbf}); // vals of Sij
    
    {
        xt::xtensor<double, 2> sij_mutable{sij}; // sij is const, make a deep copy
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

 
    // 依据初始猜测构建 Fock 矩阵
    // fprime @ cprime = cprime @ e
    // sij_inv_half @ cprime = vecs
    xt::xtensor<double, 2> fprime = xt::zeros<double>({nbf, nbf}); // fprime = (sij_inv_half @ hij) @ sij_inv_half ; will be reused in eigh
    xt::xtensor<double, 2> cprime = xt::zeros<double>({nbf, nbf}); // will be reused
    xt::xtensor<double, 1> e = xt::zeros<double>({nbf});           // will be reused
    xt::xtensor<double, 2> vecs = xt::zeros<double>({nbf, nbf});   // will be reused
    if (initial_guess == "core") {
        xt::xtensor<double, 2> tmp = xt::zeros<double>({nbf, nbf});
        cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                    CblasNoTrans, nbf, nbf, 
                    nbf, 1., sij_inv_half.data(), 
                    nbf, hij.data(), nbf, 
                    0., tmp.data(), nbf); // sij_inv_half @ hij = tmp
        cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                    CblasNoTrans, nbf, nbf, 
                    nbf, 1., tmp.data(), 
                    nbf, sij_inv_half.data(), nbf, 
                    0., fprime.data(), nbf); // fprime = (sij_inv_half @ hij) @ sij_inv_half

        xt::xtensor<double, 2> fprime_mutable{fprime}; // let fprime_mutable be destroyed
        double sfmin = LAPACKE_dlamch('S'); // safe minium for DSYEVR
        int m;
        std::vector<int> isuppz(2*nbf, 0);
        int info = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, 'V', 'I', 'U', 
                                nbf, fprime_mutable.data(), nbf, 0., 
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
    else if(initial_guess == "SAD") {std::cerr << "SAD not implemented yet\n"; exit(-1);}
    else {std::cerr << "Unknown initial guess\n"; exit(-1);}

    xt::xtensor<double, 2> Puv = xt::zeros<double>({nbf, nbf}); // Puv = vecs[:,:nocc] @ vecs[:,:nocc].T
    cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                CblasTrans, nbf, nbf, 
                this->nocc, 2., vecs.data(), 
                nbf, vecs.data(), nbf, 
                0., Puv.data(), nbf);
    


    // 初始化 libxc
    xc_func_type funcx, funcc;
    xc_func_init(&funcx, X_id, XC_UNPOLARIZED); // exchange
    xc_func_init(&funcc, C_id, XC_UNPOLARIZED); // correlation


    ///////////////////
    //>! SCF START !<//
    //////////////////
    std::cout << std::string(65, '=') << std::endl;
    std::cout << std::string(30, '=') << " SCF " << std::string(30, '=') << std::endl;
    std::cout << std::string(65, '=') << std::endl;
    std::cout << "bf type " << (pure ? "pure" : "cartesian") << std::endl;
    std::cout << "Fock shape " << "(" << hij.shape(0) << ", " << hij.shape(1) << ")" << std::endl;

    while (true) {
    std::cout << "\n!>STEP " << counter+1 << std::endl;
    xt::xtensor<double, 2> rho = xt::zeros<double>({this->natom, this->natgrid});

#ifdef TIMEIT
    const auto t0{std::chrono::steady_clock::now()};
#endif

#pragma omp parallel for
    for (auto iatom=0; iatom < this->natom; ++iatom) {
        for (auto ibf=0; ibf < nbf; ++ibf) {
            for (auto jbf=ibf; jbf < nbf; ++jbf) {
                auto phi_mu = xt::row(aos_vals[iatom], ibf);
                auto phi_nu = xt::row(aos_vals[iatom], jbf);
                if (jbf > ibf) xt::row(rho, iatom) += 2. * phi_mu * Puv(ibf, jbf) * phi_nu;
                else           xt::row(rho, iatom) +=      phi_mu * Puv(ibf, jbf) * phi_nu;
            }
        }
    }

#ifdef TIMEIT
    const auto t1{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> wtime_rho{t1 - t0};
#endif


    double ne_quad = xt::sum(rho * mweights)(); // 检查通过求积得到的电子数是否复现真值
    std::cout << std::format("Number of electron via quadrature = {:<.15f}\n", ne_quad);
    rho *= this->ne / ne_quad;

#ifdef TIMEIT
    const auto t0_Juv{std::chrono::steady_clock::now()};
#endif
    xt::xtensor<double, 2> Juv  = pure ? sf::get_J(this->mol.shells_pure, Puv) : sf::get_J(this->mol.shells_cart, Puv);
#ifdef TIMEIT
    const auto t1_Juv{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> wtime_Juv{t1_Juv - t0_Juv};
#endif

    xt::xtensor<double, 2> Kuv  = xt::zeros<double>({nbf, nbf});
    xt::xtensor<double, 2> Cuv  = xt::zeros<double>({nbf, nbf});
    xt::xtensor<double, 2> Exuv = xt::zeros<double>({nbf, nbf});
    xt::xtensor<double, 2> Ecuv = xt::zeros<double>({nbf, nbf});
    xt::xtensor<double, 2> vx   = xt::zeros<double>({this->natom, this->natgrid});
    xt::xtensor<double, 2> ex   = xt::zeros<double>({this->natom, this->natgrid});
    xt::xtensor<double, 2> vc   = xt::zeros<double>({this->natom, this->natgrid});
    xt::xtensor<double, 2> ec   = xt::zeros<double>({this->natom, this->natgrid});

#ifdef TIMEIT
    const auto t2{std::chrono::steady_clock::now()};
#endif
    xc_lda_vxc(&funcx, this->ngrid, rho.data(), vx.data()); // this is why vx cannot be vector<xtensor<double, 1>>
    xc_lda_exc(&funcx, this->ngrid, rho.data(), ex.data());
    xc_lda_vxc(&funcc, this->ngrid, rho.data(), vc.data());
    xc_lda_exc(&funcc, this->ngrid, rho.data(), ec.data());
#ifdef TIMEIT
    const auto t3{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> wtime_libxc{t3 - t2};
#endif


// rho (natom, natgrid)
// aos_vals2 [natom, (nbf, natgrid)]
// mweights (natom, natgrid)
#ifdef TIMEIT
    const auto t4{std::chrono::steady_clock::now()};
#endif


// BLAS level 3 for fock is 10x faster
#pragma omp parallel for
    for (auto iatom=0; iatom < this->natom; ++iatom) {
        xt::xtensor<double, 1> vx_mweights = xt::row(mweights, iatom) * xt::row(vx, iatom);
        xt::xtensor<double, 1> vc_mweights = xt::row(mweights, iatom) * xt::row(vc, iatom);
        xt::xtensor<double, 1> ex_mweights = xt::row(mweights, iatom) * xt::row(ex, iatom);
        xt::xtensor<double, 1> ec_mweights = xt::row(mweights, iatom) * xt::row(ec, iatom);
        xt::xtensor<double, 2> vx_weighted_aos_vals = aos_vals[iatom] * vx_mweights;
        xt::xtensor<double, 2> vc_weighted_aos_vals = aos_vals[iatom] * vc_mweights;
        xt::xtensor<double, 2> ex_weighted_aos_vals = aos_vals[iatom] * ex_mweights;
        xt::xtensor<double, 2> ec_weighted_aos_vals = aos_vals[iatom] * ec_mweights;
        xt::xtensor<double, 2> Kuv_local  = xt::zeros<double>({nbf, nbf});
        xt::xtensor<double, 2> Cuv_local  = xt::zeros<double>({nbf, nbf});
        xt::xtensor<double, 2> Exuv_local = xt::zeros<double>({nbf, nbf});
        xt::xtensor<double, 2> Ecuv_local = xt::zeros<double>({nbf, nbf});
        cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                    CblasTrans, nbf, nbf, 
                    this->natgrid, 1., vx_weighted_aos_vals.data(), 
                    this->natgrid, aos_vals[iatom].data(), this->natgrid, 
                    0., Kuv_local.data(), nbf); // (vx x mweights x aos_vals) @ aos_vals.T => Kuv

        cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                    CblasTrans, nbf, nbf, 
                    this->natgrid, 1., vc_weighted_aos_vals.data(), 
                    this->natgrid, aos_vals[iatom].data(), this->natgrid, 
                    0., Cuv_local.data(), nbf); // (vc x mweights x aos_vals) @ aos_vals.T => Cuv

        cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                    CblasTrans, nbf, nbf, 
                    this->natgrid, 1., ex_weighted_aos_vals.data(), 
                    this->natgrid, aos_vals[iatom].data(), this->natgrid, 
                    0., Exuv_local.data(), nbf); // (ex x mweights x aos_vals) @ aos_vals.T => Exuv

        cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                    CblasTrans, nbf, nbf, 
                    this->natgrid, 1., ec_weighted_aos_vals.data(), 
                    this->natgrid, aos_vals[iatom].data(), this->natgrid, 
                    0., Ecuv_local.data(), nbf); // (ec x mweights x aos_vals) @ aos_vals.T => Ecuv
#pragma omp critical
{
    Kuv  += Kuv_local;
    Cuv  += Cuv_local;
    Exuv += Exuv_local;
    Ecuv += Ecuv_local;
}
    }


#ifdef TIMEIT
    const auto t5{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> wtime_fock{t5 - t4};
#endif



    xt::xtensor<double, 2> fock = tij + vij + Juv + Kuv + Cuv;
    fock = (fock + xt::transpose(fock)) * 0.5;
    focks.emplace_back(fock);


    // 能量分解
    // rho (natom, natgrid)
    // mweights (natom, natgrid)
    ne_quad              = xt::sum(rho * mweights)();
    xt::xtensor<double, 2> Puv_times_tijT  = xt::zeros<double>({nbf, nbf});
    xt::xtensor<double, 2> Puv_times_vijT  = xt::zeros<double>({nbf, nbf});
    xt::xtensor<double, 2> Puv_times_JuvT  = xt::zeros<double>({nbf, nbf});
    xt::xtensor<double, 2> Puv_times_ExuvT = xt::zeros<double>({nbf, nbf});
    xt::xtensor<double, 2> Puv_times_EcuvT = xt::zeros<double>({nbf, nbf});

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
                nbf, Exuv.data(), nbf, 
                0., Puv_times_ExuvT.data(), nbf); // Puv @ Exuv.T
    cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                CblasTrans, nbf, nbf, 
                nbf, 1., Puv.data(), 
                nbf, Ecuv.data(), nbf, 
                0., Puv_times_EcuvT.data(), nbf); // Puv @ Ecuv.T
    

    double kin_e         = ltr(Puv_times_tijT.data(),  nbf);
    double ext_e         = ltr(Puv_times_vijT.data(),  nbf);
    double hartree_e     = ltr(Puv_times_JuvT.data(),  nbf) * 0.5;
    double exchange_e    = ltr(Puv_times_ExuvT.data(), nbf);
    double correlation_e = ltr(Puv_times_EcuvT.data(), nbf);
    double etot = kin_e + ext_e + hartree_e + exchange_e + correlation_e + this->mol.e_nuc;

    this->ne_quads.emplace_back(ne_quad);
    this->kinetic_energies.emplace_back(kin_e);
    this->external_energies.emplace_back(ext_e);
    this->hartree_energies.emplace_back(hartree_e);
    this->exchange_energies.emplace_back(exchange_e);
    this->correlation_energies.emplace_back(correlation_e);
    this->etots.emplace_back(etot);
    this->energy_log();



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

#ifdef TIMEIT
    const auto t6{std::chrono::steady_clock::now()};
#endif
    
    if (counter >= 2) cDIIS(fock, focks, diis_res);

#ifdef TIMEIT
    const auto t7{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> wtime_diis{t7 - t6};
#endif
    
#ifdef TIMEIT
    const auto t8{std::chrono::steady_clock::now()};
#endif
    // 对角化 Fock matrix
    {
        xt::xtensor<double, 2> tmp = xt::zeros<double>({nbf, nbf});
        cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                    CblasNoTrans, nbf, nbf, 
                    nbf, 1., sij_inv_half.data(), 
                    nbf, fock.data(), nbf, 
                    0., tmp.data(), nbf); // sij_inv_half @ fock = tmp
        cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                    CblasNoTrans, nbf, nbf, 
                    nbf, 1., tmp.data(), 
                    nbf, sij_inv_half.data(), nbf, 
                    0., fprime.data(), nbf); // (sij_inv_half @ fock) @ sij_inv_half = fprime

        xt::xtensor<double, 2> fprime_mutable{fprime}; // let fprime_mutable be destroyed
        double sfmin = LAPACKE_dlamch('S'); // safe minium for DSYEVR
        int m;
        std::vector<int> isuppz(2*nbf, 0);
        int info = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, 'V', 'I', 'U', 
                                nbf, fprime_mutable.data(), nbf, 0., 
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

    cblas_dgemm(CblasRowMajor, CblasNoTrans, 
                CblasTrans, nbf, nbf, 
                this->nocc, 2., vecs.data(), 
                nbf, vecs.data(), nbf, 
                0., Puv.data(), nbf); // Puv = 2 x vecs[:,:nocc] @ vecs[:,:nocc].T

#ifdef TIMEIT
    const auto t9{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> wtime_diag{t9 - t8};
#endif


    // 收敛判断
    double e_diff = std::abs(this->etots.back() - this->etots[this->etots.size()-2]);
    double d_diff = xt::sum(xt::square(Puv - Puv_old))();
    if (e_diff < e_convergence && d_diff < d_convergence) {
        std::cout << std::format("Density difference = {:<.12f}\n", d_diff);
        std::cout << std::format("SCF converged after {:d} step. Total energy = {:<.15f}\n", counter+1, this->etots.back());
        
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
        std::cout << "\n!> NORMAL TERMINATION" << std::endl;
    
        break;
    }

    if (counter >= maxiter) {
        std::cout << "SCF did not converge in " << maxiter << " step. Energy in the last iteration " << this->etots.back() << std::endl;
        break;
    }

    std::cout << std::format("Density difference = {:.12f}\n", d_diff);

    counter++;
    Puv_old = Puv;

#ifdef TIMEIT
    const auto t_while_end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> wtime_while{t_while_end - t0};
    std::cout << "wtime Juv   : " << wtime_Juv   << std::endl;
    std::cout << "wtime rho   : " << wtime_rho   << std::endl;
    std::cout << "wtime libxc : " << wtime_libxc << std::endl;
    std::cout << "wtime fock  : " << wtime_fock  << std::endl;
    // std::cout << "wtime diis  : " << wtime_diis << std::endl;
    std::cout << "wtime diag  : " << wtime_fock  << std::endl;
    std::cout << "wtime       : " << wtime_while << std::endl;
#endif

    } // end scf while


    xc_func_end(&funcx);
    xc_func_end(&funcc);

    std::cout << "safe here\n";
} // end DFT.scf()

