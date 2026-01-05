#include "../include/dft.hpp"
#include "xtensor-blas/xlinalg.hpp" // icpx -std=c++20 src.cc -DHAVE_CBLAS=1 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm -ldl for optimal xlinalg performance

#ifndef RELEASE
#define RELEASE
#endif


sf::DFT::DFT::DFT(const std::string& xyzfile, const std::string& name) : xyzfile(xyzfile), name(name)
{
    this->mol = sf::Molecule{xyzfile, name};
}

void sf::DFT::DFT::header_log() const
{

    // 在开始 SCF 前输出一些信息
    std::cout << "Number of atoms          " << this->natom << std::endl;
    std::cout << "Number of electrons      " << this->ne    << std::endl;
    // std::cout << "Number of aos            " << this->nao   << std::endl;
    std::cout << "Number of radial grid    " << this->nrad  << std::endl;
    std::cout << "Number of angular grid   " << this->nang  << std::endl;
    std::cout << "Number of total grid     " << this->ngrid << std::endl;
    std::cout << "Maximum angular momentum " << this->lmax  << std::endl;
    std::cout << "Nuclear repulsion energy " << std::setprecision(15) << this->mol.e_nuc << " a.u." << std::endl;
    std::cout << std::endl;
}

void sf::DFT::DFT::energy_log() const
{
    // 输出 SCF 每一步的能量分解
    std::cout << "Number of electron via quadrature " << std::fixed << std::setprecision(12) << std::left << this->ne_quads.back() << std::endl;
    std::cout << "One-electron energy               " << std::fixed << std::setprecision(12) << std::left << this->kinetic_energies.back() + this->external_energies.back() << std::endl;
    std::cout << "Two-electron energy               " << std::fixed << std::setprecision(12) << std::left << this->hartree_energies.back() << std::endl;
    std::cout << "Exchange-correlation energy       " << std::fixed << std::setprecision(12) << std::left << this->exchange_energies.back() + this->correlation_energies.back() << "    " << this->exchange_energies.back() << "    " << this->correlation_energies.back() << std::endl;
    std::cout << "Total energy                      " << std::fixed << std::setprecision(12) << std::left << this->etots.back() << std::endl;
    std::cout << "Energy difference                 " << std::fixed << std::setprecision(12) << std::left << this->etots.back() - this->etots[this->etots.size()-2] << std::endl;
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
    this->lmax     = std::floor(angular_level / 2.); // for spherical expansion truncation, not for integral evaluation!
    this->nlm      = (this->lmax + 1) * (this->lmax + 1); // l=0 -> 1 ; l=2 -> 9 ; l=14 -> 225 etc.
    this->lm.reserve(this->nlm); // (l, m) pairs
    for (auto l=0; l <= this->lmax; ++l) {
        for (auto m=-l; m <= l; ++m) {
            this->lm.emplace_back(l, m);
        }
    }

    auto it = std::ranges::find(LEBEDEV_ORDER, angular_level);
    assert(it != std::end(LEBEDEV_ORDER));
    int index = it - std::begin(LEBEDEV_ORDER);

    this->nang    = LEBEDEV_LEVEL[index];                  // 角度网格数
    this->nrad    = static_cast<int>(becke.ncheb);         // 径向网格数
    this->natgrid = this->nrad * this->nang;               // 原子网格数
    this->ngrid   = this->natom * this->nrad * this->nang; // 总网格数

    
    // grid (natom, nrad, nang, 3+natom) 坐标+权重
    // grid_global (natom, natgrid, 3) 仅坐标
    // grid_total (ngrid, 3) 全部的坐标
    // `BeckeFuzzyCell::build_grid` 未调用的话后续的初始化无法进行
    // 构建 this->aos_vals_cart {nbf_cart, natom, natgrid}
    const auto [grid, grid_global] = becke.build_grid();
    const auto grid_total = xt::reshape_view(grid_global, {ngrid, 3});
    const auto& grid_global_ref = grid_global;
    this->aos_vals_cart.resize({static_cast<size_t>(this->nbf_cart), static_cast<size_t>(this->natom), static_cast<size_t>(this->natgrid)});
    this->aos_vals_pure.resize({static_cast<size_t>(this->nbf_pure), static_cast<size_t>(this->natom), static_cast<size_t>(this->natgrid)});
    
    sf::BFs bfs{this->mol};
    const xt::xtensor<double, 2> TF = this->mol.make_tf_xt();

    #pragma omp parallel for schedule(dynamic)
    for (auto iatom=0; iatom < this->natom; ++iatom) {
        auto x_view = xt::view(grid_global_ref, iatom, xt::all(), 0); // x
        auto y_view = xt::view(grid_global_ref, iatom, xt::all(), 1); // y
        auto z_view = xt::view(grid_global_ref, iatom, xt::all(), 2); // z
        auto view_cart = xt::view(this->aos_vals_cart, xt::all(), iatom, xt::all());
        view_cart = bfs.ao_val(x_view, y_view, z_view); // {nbf_cart, natgrid}
        
        // 构建球谐形式的基函数在网格上的值矩阵
        // aos_vals_cart {nbf_cart, natom, natgrid} -> aos_vals_pure {nbf_pure, natom, natgrid}
        xt::view(this->aos_vals_pure, xt::all(), iatom, xt::all()) = xt::linalg::dot(TF, xt::xtensor<double, 2>{view_cart});
    }
    
    
    // 构建全局距离和固体角
    // dist, theta, phi (natom, ngrid)
    this->dist.resize({static_cast<size_t>(this->natom),  static_cast<size_t>(this->ngrid)}); // (natom. ngrid)
    this->theta.resize({static_cast<size_t>(this->natom), static_cast<size_t>(this->ngrid)}); // (natom. ngrid)
    this->phi.resize({static_cast<size_t>(this->natom),   static_cast<size_t>(this->ngrid)}); // (natom. ngrid)
    for (auto iatom=0; iatom < this->natom; ++iatom) {
        xt::xtensor<double, 2> _disp = grid_total - xt::row(this->mol.xyz, iatom);
        xt::row(dist, iatom)     = xt::norm_l2(_disp, {1});
        _disp /= xt::view(xt::row(this->dist, iatom), xt::all(), xt::newaxis()); // normalized
        xt::row(this->theta, iatom) = xt::acos(xt::col(_disp, 2));
        xt::row(this->phi,   iatom) = xt::atan2(xt::col(_disp, 1), xt::col(_disp, 0));
    }

    // 构建实球谐函数 ylm (natom, nlm, ngrid)
    this->ylm.resize({static_cast<size_t>(this->natom), static_cast<size_t>(this->nlm), static_cast<size_t>(this->ngrid)});
#pragma omp parallel for schedule(dynamic) collapse(2)
    for (auto iatom=0; iatom < this->natom; ++iatom) {
        for (auto ilm=0; ilm < this->nlm; ++ilm) {
            for (auto igrid=0; igrid < this->ngrid; ++igrid) {
                this->ylm(iatom, ilm, igrid) = rsh::rsh(this->lm[ilm].first, this->lm[ilm].second, this->theta(iatom, igrid), this->phi(iatom, igrid));
            }
        }
    }

    // 构建单位球面上的球谐基，所有原子共用一套
    // y_jk (nang, nlm)
    xt::xtensor<double, 1> _theta = xt::acos(xt::row(this->becke.xwleb, 2));
    xt::xtensor<double, 1> _phi   = xt::atan2(xt::row(this->becke.xwleb, 1), xt::row(this->becke.xwleb, 0));
    this->y_jk.resize({static_cast<size_t>(this->nang), static_cast<size_t>(this->nlm)});
    for (auto ilm=0; ilm < this->nlm; ++ilm) {
        for (auto iang=0; iang < this->nang; ++iang) {
            y_jk(iang, ilm) = rsh::rsh(this->lm[ilm].first, this->lm[ilm].second, _theta[iang], _phi[iang]);
        }
    }
} // end DFT.init()

xt::xtensor<double, 1> sf::DFT::DFT::build_hartree_potential(const xt::xtensor<double, 2>& rho, const xt::xtensor<double, 2>& mweights, const xt::xtensor<double, 2>& zz)
{
    // rho -> rho_ik (natom, nrad, nlm)
    xt::xtensor<double, 3> rho_ik      = xt::zeros<double>({this->natom, this->nrad, this->nlm});
    // rho_ik -> u_ik (natom, nrad, nlm)
    xt::xtensor<double, 3> u_ik        = xt::zeros<double>({this->natom, this->nrad, this->nlm});
    // u_ik -> u_ik_interp (natom, nlm, ngrid)
    xt::xtensor<double, 3> u_ik_interp = xt::zeros<double>({this->natom, this->nlm, this->ngrid});
    // u_ik_interp -> u_ij (ngrid, ) ; TO BE RETURNED
    xt::xtensor<double, 1> u_ij        = xt::zeros<double>({this->ngrid});

    for (auto iatom=0; iatom < this->natom; ++iatom) {
        // 构建电子密度球谐展开系数 rho -> rho_ik
        // lhs (nrad, nang) ; rhs (nang, nlm) -> (nrad, nlm)
        xt::xtensor<double, 2> _lhs = xt::reshape_view(xt::row(this->becke.weights, iatom) * xt::row(rho, iatom), {this->nrad, this->nang});
        xt::xtensor<double, 2> _rhs = this->y_jk * xt::view(this->becke.xwleb, 3, xt::all(), xt::newaxis());
        Eigen::Map<xxd> _lhs_map{_lhs.data(), this->nrad, this->nang};
        Eigen::Map<xxd> _rhs_map{_rhs.data(), this->nang, this->nlm};
        xxd            _rho_ik = _lhs_map * _rhs_map; // (nrad, nlm)
        xt::view(rho_ik, iatom, xt::all(), xt::all()) = xt::adapt(_rho_ik.data(), _rho_ik.size(), xt::no_ownership(), xt::xtensor<double, 2>::shape_type{static_cast<size_t>(this->nrad), static_cast<size_t>(this->nlm)});

        // 求解泊松方程 rho_ik -> u_ik
        double qn = xt::sum(xt::row(mweights, iatom) * xt::row(rho, iatom))(); // partial charge
#pragma omp parallel for schedule(dynamic)
        for (auto ilm=0; ilm < this->nlm; ++ilm) {
            auto _rho_ik_view           = xt::view(rho_ik, iatom, xt::all(), ilm); // (nrad, )
            auto _u_ik_view             = xt::view(u_ik,   iatom, xt::all(), ilm); // (nrad, )
            xt::xtensor<double, 1> _poisson = beckegrid::poisson_solver(_rho_ik_view, this->lm[ilm].first, xt::row(this->becke.rcheb, iatom), this->mol.rms[iatom], qn);
            _u_ik_view = xt::view(_poisson, xt::range(1, this->nrad+1)); // [1:-1] (nrad, )

            // 给定一个 lm, 插值计算所有网格处的库伦势球谐展开系数 u_ik_interp
            // cspline::CubicSpline cs(this->becke.zcheb, _u_ik_view);
            // cs.generate_spline();
            gsl_interp_accel*     acc = gsl_interp_accel_alloc();
            gsl_spline*        spline = gsl_spline_alloc(gsl_interp_cspline, this->nrad); // bc = 'natural'
            xt::xtensor<double, 1> _y = xt::view(_poisson, xt::range(1, this->nrad+1));
            gsl_spline_init(spline, this->becke.zcheb.data(), _y.data(), this->nrad);
                for (auto igrid=0; igrid < this->ngrid; ++igrid) { // 不要对这层循环使用 omp
                    if (zz(iatom, igrid) > this->becke.zcheb.back()) {
                        u_ik_interp(iatom, ilm, igrid) = 0.;   // z = nrad <=> r=0
                    }
                    else if (zz(iatom, igrid) < this->becke.zcheb.front()) { // z = 1 <=> r->inf
                        if (ilm == 0) u_ik_interp(iatom, ilm, igrid) = std::sqrt(4. * PI) * qn;
                        else          u_ik_interp(iatom, ilm, igrid) = 0.;
                    }
                    // else              u_ik_interp(iatom, ilm, igrid) = cs.eval(zz(iatom, igrid));
                    else              u_ik_interp(iatom, ilm, igrid) = gsl_spline_eval(spline, zz(iatom, igrid), acc);
                } // loop ngrid
            gsl_spline_free(spline);
            gsl_interp_accel_free(acc);
        } // loop lm
        u_ij += xt::sum(xt::view(u_ik_interp, iatom, xt::all(), xt::all()) / xt::row(this->dist, iatom) * xt::view(this->ylm, iatom, xt::all(), xt::all()), {0}); // (ngrid, )
    } // loop atoms

    return u_ij;
} // end buid_hartree_potential

void sf::DFT::DFT::diis(xxd& fock, const std::vector<xxd>& focks, const std::vector<xxd>& diis_res)
{
    std::cout << "DIIS enabled\n";
    xxd B = xxd::Zero(focks.size()+1, focks.size()+1);
    B.bottomRows(1)         = -Eigen::RowVectorXd::Ones(focks.size()+1); // B[-1, :] = -1
    B.rightCols(1)          = -xd::Ones(focks.size()+1); // B[:, -1] = -1
    B(B.rows()-1, B.cols()-1) = 0.; // B[-1, -1] = 0.

    for (auto ii=0; ii < focks.size(); ++ii) {
        for (auto jj=0; jj < focks.size(); ++jj) {
            B(ii, jj) = (diis_res[ii] * diis_res[jj].transpose()).trace();
        }
    }

    xd diis_rhs = xd::Zero(B.rows());
    diis_rhs[diis_rhs.size()-1] = -1.; // diis_res[-1] = -1.
    xd diis_coef = B.lu().solve(diis_rhs);
    fock *= 0.;
    for (auto iii=0; iii < diis_coef.size()-1; ++iii) {
        fock += focks[iii] * diis_coef[iii];
    }
    fock = (fock + fock.transpose()) * 0.5;
}

void sf::DFT::DFT::scf(const int maxiter, 
                       const double e_convergence, 
                       const double d_convergence, 
                       const int nbuffer, 
                       const std::string &initial_guess, 
                       const int X_id, 
                       const int C_id, 
                       const bool pure)
{
    // 单电子算符 (nao, nao)
    const std::vector<libint2::Atom>& atoms = this->mol.atoms;
    const xt::xtensor<double, 3>& aos_vals = pure ? this->aos_vals_pure : this->aos_vals_cart; // (nbf, natom, natgrid)
    int nbf = pure ? this->nbf_pure : this->nbf_cart;

    const xt::xtensor<double, 2> sij = pure ? sf::get_olp(this->mol.shells_pure)        : sf::get_olp(this->mol.shells_cart);
    const xt::xtensor<double, 2> tij = pure ? sf::get_kin(this->mol.shells_pure)        : sf::get_kin(this->mol.shells_cart);
    const xt::xtensor<double, 2> vij = pure ? sf::get_ext(this->mol.shells_pure, atoms) : sf::get_ext(this->mol.shells_cart, atoms);
    assert(sij.shape(0) == nbf && sij.shape(1) == nbf);
    const Eigen::Map<const xxd> sij_map{sij.data(), nbf, nbf};
    const Eigen::Map<const xxd> tij_map{tij.data(), nbf, nbf};
    const Eigen::Map<const xxd> vij_map{vij.data(), nbf, nbf};
    const xxd hij = tij_map + vij_map;

   
    // mweights (natom, natgrid)
    // 径向权重（r^2 wrad），角度权重（wleb），网格权重（becke.weights (natom, natgrid)）
    xt::xtensor<double, 2> mweights = xt::zeros<double>({this->natom, this->natgrid});
    for (auto iatom=0; iatom < this->natom; ++iatom) {
        xt::xtensor<double, 1> _wrad = xt::row(this->becke.rcheb, iatom) * xt::row(this->becke.rcheb, iatom) * xt::row(this->becke.wcheb, iatom); // (nrad, ) 径向权重
        xt::xtensor<double, 1> _wang = xt::row(this->becke.xwleb, 3); // (nang, ) 角度权重
        Eigen::Map<xd> _wrad_map{_wrad.data(), static_cast<int>(_wrad.size())}; // 用来做外积
        Eigen::Map<xd> _wang_map{_wang.data(), static_cast<int>(_wang.size())}; // 用来做外积
        xxd _outer = _wrad_map * _wang_map.transpose();
        auto _outer_adaptor = xt::adapt(_outer.data(), _outer.size(), xt::no_ownership(), xt::xtensor<double, 2>::shape_type{static_cast<size_t>(this->nrad), static_cast<size_t>(this->nang)});
        xt::row(mweights, iatom) = xt::ravel(_outer_adaptor) * xt::row(this->becke.weights, iatom);
    }

    // zz (natom, ngrid)
    // z -> x -> r, u_ik = u_ik(r[x[z]]) 通过内插 z 计算全局 u_ik，进而构建全局库伦势 u_ij
    xt::xtensor<double, 2> zz = xt::zeros<double>({this->natom, this->ngrid});
    const double _prefactor = (this->nrad + 1) / PI;
    for (auto iatom=0; iatom < this->natom; ++iatom) {
        auto _dist = xt::row(this->dist, iatom);
        xt::row(zz, iatom) = _prefactor * xt::acos((_dist - this->mol.rms[iatom]) / (_dist + this->mol.rms[iatom]));
    }

    // 输出一些信息
    this->header_log();

    // 参数设置
    int counter = 0;
    this->etots.emplace_back(std::nan("1"));
    xxd Puv_old = xxd::Ones(nbf, nbf) * std::nan("1");
    std::vector<xxd> focks   ; focks.reserve(maxiter + 1);
    std::vector<xxd> diis_res; diis_res.reserve(maxiter + 1);

    // 重叠矩阵正交化
    Eigen::SelfAdjointEigenSolver<xxd> eigsolver; // 厄米矩阵特征值分解
    eigsolver.compute(sij_map);
    assert(eigsolver.info() == Eigen::Success);
    const xxd u = eigsolver.eigenvectors();
    const xd  s = eigsolver.eigenvalues();
    const xxd inv_s = (1. / s.array()).sqrt().matrix().asDiagonal().toDenseMatrix(); // s^-1/2
    const xxd sij_inv_half = u * inv_s * u.transpose();

    // 依据初始猜测构建 Fock 矩阵
    xxd fprime = xxd::Zero(nbf, nbf);
    xxd cprime = xxd::Zero(nbf, nbf);
    xd  e      =  xd::Zero(nbf);
    xxd vecs   = xxd::Zero(nbf, nbf);
    if (initial_guess == "core") {
        fprime = sij_inv_half * hij * sij_inv_half; // hij == fock
        eigsolver.compute(fprime);
        assert(eigsolver.info() == Eigen::Success);
        cprime = eigsolver.eigenvectors(); // C'
        e      = eigsolver.eigenvalues();  // 升序排列？
        vecs   = sij_inv_half * cprime;    // C
    }
    else {std::cerr << "not implemented yet\n"; exit(-1);}

    // 构建初始密度矩阵
    // Puv = C[:, :nocc] @ C[:, :nocc].T
    auto block = vecs.block(0, 0, nbf, this->nocc);
    xxd Puv = 2. * block * block.transpose();

    //>!
    std::cout << std::string(65, '=') << std::endl;
    std::cout << std::string(30, '=') << " SCF " << std::string(30, '=') << std::endl;
    std::cout << std::string(65, '=') << std::endl;
    //>!

    this->kinetic_energies.reserve(maxiter + 1);
    this->external_energies.reserve(maxiter + 1);
    this->hartree_energies.reserve(maxiter + 1);
    this->exchange_energies.reserve(maxiter + 1);
    this->correlation_energies.reserve(maxiter + 1);
    this->etots.reserve(maxiter + 1);
    this->ne_quads.reserve(maxiter + 1);

    // 初始化 libxc
    xc_func_type funcx, funcc;
    xc_func_init(&funcx, X_id, XC_UNPOLARIZED); // exchange
    xc_func_init(&funcc, C_id, XC_UNPOLARIZED); // correlation


    ///////////////////
    //>! SCF START !<//
    //////////////////
    std::cout << "bf type " << (pure ? "pure" : "cartesian") << std::endl;
    std::cout << "Fock shape " << "(" << hij.rows() << ", " << hij.cols() << ")" << std::endl;
#ifdef RELEASE
    // DEBUG 模式 SCF 只运行一步
    while (true) {
#endif

    std::cout << "!>STEP " << counter+1 << std::endl;
    // 更新电子密度 rho (natom, natgrid)
    // 构建电子密度 Puv -> rho use aos_vals_pure if PURE
    xt::xtensor<double, 2> rho = xt::zeros<double>({this->natom, this->natgrid});
#pragma omp parallel for
    for (auto iatom=0; iatom < this->natom; ++iatom) {
        for (auto ibf=0; ibf < nbf; ++ibf) {
            for (auto jbf=ibf; jbf < nbf; ++jbf) {
                auto _phi_mu = xt::view(aos_vals, ibf, iatom, xt::all());
                auto _phi_nu = xt::view(aos_vals, jbf, iatom, xt::all());
                if (jbf > ibf) xt::row(rho, iatom) += 2. * _phi_mu * Puv(ibf, jbf) * _phi_nu;
                else           xt::row(rho, iatom) +=      _phi_mu * Puv(ibf, jbf) * _phi_nu;
            }
        }
    }

    // 检查通过求积得到的电子数是否复现真值
    double ne_quad = xt::sum(rho * mweights)();
    std::cout << "Number of electron via quadrature " << std::setprecision(15) << std::fixed << ne_quad << std::endl;
    // assert(std::abs(ne_quad - this->ne) < 1e-3);
    rho *= this->ne / ne_quad;

    // 构建全局库伦势
    auto u_ij= build_hartree_potential(rho, mweights, zz);

    // 构建 Juv, Kuv, Cuv, Exuv, Ecuv ; 调用 libxc
    xxd Juv  = xxd::Zero(nbf, nbf);
    xxd Kuv  = xxd::Zero(nbf, nbf);
    xxd Cuv  = xxd::Zero(nbf, nbf);
    xxd Exuv = xxd::Zero(nbf, nbf);
    xxd Ecuv = xxd::Zero(nbf, nbf);

    auto rho_flattened        = xt::ravel(rho); // (natom, natgrid) -> (ngrid, )
    xt::xtensor<double, 1> vx = xt::zeros<double>({this->ngrid});
    xt::xtensor<double, 1> ex = xt::zeros<double>({this->ngrid});
    xt::xtensor<double, 1> vc = xt::zeros<double>({this->ngrid});
    xt::xtensor<double, 1> ec = xt::zeros<double>({this->ngrid});

    xc_lda_vxc(&funcx, this->ngrid, rho_flattened.data(), vx.data());
    xc_lda_exc(&funcx, this->ngrid, rho_flattened.data(), ex.data());
    xc_lda_vxc(&funcc, this->ngrid, rho_flattened.data(), vc.data());
    xc_lda_exc(&funcc, this->ngrid, rho_flattened.data(), ec.data());

#pragma omp parallel for
    for (auto ibf=0; ibf < nbf; ++ibf) {
        for (auto jbf=ibf; jbf < nbf; ++jbf) {
            auto phi_mu_flattened = xt::ravel(xt::view(aos_vals, ibf, xt::all(), xt::all())); // (ngrid, )
            auto phi_nu_flattened = xt::ravel(xt::view(aos_vals, jbf, xt::all(), xt::all())); // (ngrid, )
            auto mweights_flattened = xt::ravel(mweights); // (ngrid, )
            Juv(ibf,  jbf) = xt::sum(phi_mu_flattened * u_ij * phi_nu_flattened * mweights_flattened)();
            Kuv(ibf,  jbf) = xt::sum(phi_mu_flattened * vx   * phi_nu_flattened * mweights_flattened)();
            Cuv(ibf,  jbf) = xt::sum(phi_mu_flattened * vc   * phi_nu_flattened * mweights_flattened)();
            Exuv(ibf, jbf) = xt::sum(phi_mu_flattened * ex   * phi_nu_flattened * mweights_flattened)();
            Ecuv(ibf, jbf) = xt::sum(phi_mu_flattened * ec   * phi_nu_flattened * mweights_flattened)();
            if (jbf > ibf) {
                Juv(jbf,  ibf) = Juv(ibf,  jbf);
                Kuv(jbf,  ibf) = Kuv(ibf,  jbf);
                Cuv(jbf,  ibf) = Cuv(ibf,  jbf);
                Exuv(jbf, ibf) = Exuv(ibf, jbf);
                Ecuv(jbf, ibf) = Ecuv(ibf, jbf);
            }
        }
    }

    xxd fock = tij_map + vij_map + Juv + Kuv + Cuv;
    fock = (fock + fock.transpose()) * 0.5;
    focks.emplace_back(fock);
    diis_res.emplace_back(sij_inv_half * (fock * Puv * sij_map - sij_map * Puv * fock) * sij_inv_half);

    // 能量分解
    ne_quad              = xt::sum(rho * mweights)();
    double kin_e         = (Puv *  tij_map.transpose()).trace();
    double ext_e         = (Puv *  vij_map.transpose()).trace();
    double hartree_e     = (Puv *  Juv.transpose()).trace() * 0.5;
    double exchange_e    = (Puv * Exuv.transpose()).trace();
    double correlation_e = (Puv * Ecuv.transpose()).trace();
    double etot          = kin_e + ext_e + hartree_e + exchange_e + correlation_e + this->mol.e_nuc;
    this->ne_quads.emplace_back(ne_quad);
    this->kinetic_energies.emplace_back(kin_e);
    this->external_energies.emplace_back(ext_e);
    this->hartree_energies.emplace_back(hartree_e);
    this->exchange_energies.emplace_back(exchange_e);
    this->correlation_energies.emplace_back(correlation_e);
    this->etots.emplace_back(etot);
    this->energy_log();

    // 收敛加速
    if (counter > 2) diis(fock, focks, diis_res);

    // 对角化 Fock matrix
    fprime = sij_inv_half * fock * sij_inv_half;
    eigsolver.compute(fprime);
    assert(eigsolver.info() == Eigen::Success);
    cprime = eigsolver.eigenvectors();
    e      = eigsolver.eigenvalues();
    vecs   = sij_inv_half * cprime;
    // block  = vecs.block(0, 0, this->mol.nao, nocc);
    block  = vecs.block(0, 0, nbf, nocc);
    Puv    = 2. * block * block.transpose();

    // 收敛判断
    double e_diff = std::abs(this->etots.back() - this->etots[this->etots.size()-2]);
    double d_diff = std::sqrt((Puv.array() - Puv_old.array()).square().sum());
    if (e_diff < e_convergence && d_diff < d_convergence) {
        std::cout << "Density difference " << std::setprecision(12) << d_diff << std::endl;
        std::cout << std::endl;
        std::cout << "SCF converged after " << counter+1 << " step. Total energy " << std::setprecision(15) << this->etots.back() << std::endl;
        
        int iorb = 1;
        std::cout << "Doubly occupied:\n";
        for (auto val : e.head(this->nocc).transpose()) {
            std::cout << std::setw(15) << std::setprecision(10) << std::fixed << std::right << val << "    ";
            if (iorb % 4 == 0) std::cout << std::endl;
            iorb += 1;
        }
        
        iorb = 1;
        std::cout << "\nVirtual:\n";
        for (auto val : e.tail(nbf-this->nocc).transpose()) {
            std::cout << std::setw(15) << std::setprecision(10) << std::fixed << std::right << val << "    ";
            if (iorb % 4 == 0) std::cout << std::endl;
            iorb += 1;
        }
        std::cout << "\n!> NORMAL TERMINATION" << std::endl;
    
#ifdef RELEASE
        break;
#endif
    }

    if (counter >= maxiter) {
        std::cout << "SCF did not converge in " << maxiter << " step. Energy in the last iteration " << this->etots.back() << std::endl;
#ifdef RELEASE
        break;
#endif
    }

    std::cout << "Density difference        " << std::setprecision(12) << d_diff << std::endl;
    std::cout << std::endl;

    counter++;
    Puv_old = Puv;
#ifdef RELEASE
    } // end scf while
#endif

    xc_func_end(&funcx);
    xc_func_end(&funcc);

    std::cout << "safe here\n";
} // end DFT.scf()

