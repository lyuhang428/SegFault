#include "../include/hf.hpp"

#ifndef RELEASE
#define RELEASE
#endif


void sf::HF::diis(xxd& fock, const std::vector<xxd>& focks, const std::vector<xxd>& diis_res)
{
    std::cout << "    DIIS ENABLED\n";
    xxd B = xxd::Zero(focks.size()+1, focks.size()+1);
    B.bottomRows(1) = -Eigen::RowVectorXd::Ones(focks.size()+1); // B[-1, :] = -1
    B.rightCols(1)  = -xd::Ones(focks.size()+1); // B[:, -1] = -1
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




sf::HF::HF::HF(const std::string& xyzfile, const std::string& name) : xyzfile(xyzfile), name(name)
{
    // 在方法类内部初始化一个分子对象，而非传入一个分子对象引用
    this->mol = sf::Molecule{xyzfile, name};
}

void sf::HF::HF::energy_log() const
{
    // 输出 SCF 每一步的能量分解
    std::cout << "One-electron energy               " << std::fixed << std::setprecision(12) << std::left << this->kinetic_energies.back() + this->external_energies.back() << std::endl;
    std::cout << "Two-electron energy               " << std::fixed << std::setprecision(12) << std::left << this->hartree_energies.back() << std::endl;
    std::cout << "Total energy                      " << std::fixed << std::setprecision(12) << std::left << this->etots.back() << std::endl;
    std::cout << "Energy difference                 " << std::fixed << std::setprecision(12) << std::left << this->etots.back() - this->etots[this->etots.size()-2] << std::endl;
}

void sf::HF::HF::scf(int maxiter, double e_convergence, double d_convergence, int nbuffer, const std::string& initial_guess, const bool pure)
{
    // 构造函数中通过坐标文件和基组文件初始化分子对象
    // 取得分子对象中的原子和壳层信息
    // 分子类存储了纯球谐和笛卡尔两种类型的壳层，这里根据 pure 参数选择使用哪一种
    const std::vector<libint2::Atom>& atoms = this->mol.atoms;
    const std::vector<libint2::Shell>& shells = pure ? this->mol.shells_pure : this->mol.shells_cart;
    const int nbf = pure ? this->mol.nbf_pure : this->mol.nbf_cart;

    // const std::vector<libint2::Shell>& shells_pure = this->mol.shells_pure;
    // const std::vector<libint2::Shell>& shells_cart = this->mol.shells_cart;
    
    const int nshell = this->mol.nshell;
    // const int nbf_pure = this->mol.nbf_pure;
    // const int nbf_cart = this->mol.nbf_cart;
    const double e_nuc = this->mol.e_nuc;

    std::vector<size_t> shell2bf = sf::get_shell2bf(shells);  // index
    const size_t       max_nprim = sf::get_max_nprim(shells); // maximum number of primitive functions in all shell
    const int               lmax = sf::get_lmax(shells);

    // very cheap
    const xt::xtensor<double, 2> sij = sf::get_olp(shells);
    const xt::xtensor<double, 2> tij = sf::get_kin(shells);
    const xt::xtensor<double, 2> vij = sf::get_ext(shells, atoms);

    const Eigen::Map<const xxd> sij_map{sij.data(), nbf, nbf};
    const Eigen::Map<const xxd> tij_map{tij.data(), nbf, nbf};
    const Eigen::Map<const xxd> vij_map{vij.data(), nbf, nbf};
    const xxd Hcore = tij_map + vij_map;

    Eigen::SelfAdjointEigenSolver<xxd> eigsolver;
    eigsolver.compute(sij_map);
    const xxd u = eigsolver.eigenvectors();
    const xd  s = eigsolver.eigenvalues();
    const xxd inv_s = (1. / s.array()).sqrt().matrix().asDiagonal().toDenseMatrix(); // s^-1/2
    const xxd sij_inv_half = u * inv_s * u.transpose(); // for symmetrical orthogonalization and diis residue vector

    std::vector<xxd> focks   ; focks.reserve(maxiter+1);    // for diis
    std::vector<xxd> diis_res; diis_res.reserve(maxiter+1); // for diis

    // params setup
    int counter = 0;
    double etot_old = std::nan("1");
    xxd Puv_old = xxd::Ones(nbf, nbf) * std::nan("1");
    const int ne = this->mol.ne;
    const int nocc = this->mol.nocc;

    // initial guess
    // TODO: add SAD SAP
    xxd fock = Hcore;


    /////////////
    //>! SCF <!//
    /////////////
    std::cout << "bf type " << (pure ? "pure" : "cartesian") << std::endl;
    std::cout << "Fock shape " << "(" << fock.rows() << ", " << fock.cols() << ")" << std::endl;
#ifdef RELEASE
    while (true) {
#endif

        std::cout << "\n>!STEP " << counter+1 << std::endl;

        Eigen::GeneralizedSelfAdjointEigenSolver<xxd> eigsolver2;
        eigsolver2.compute(fock, sij_map);
        xd  vals = eigsolver2.eigenvalues();
        xxd vecs = eigsolver2.eigenvectors();
        xxd Puv  = 2. * vecs.leftCols(nocc) * vecs.leftCols(nocc).transpose(); // 初始密度矩阵

        // 计算库伦矩阵和交换矩阵
        // DO NOT store all ERI
        // compute (ij|kl) on the fly
        // TODO: add Schwartz screening / density fitting
        auto [_Juv, _Kuv] = sf::get_JK(shells, Puv); // 已乘密度矩阵
        Eigen::Map<xxd> Juv{_Juv.data(), nbf, nbf};
        Eigen::Map<xxd> Kuv{_Kuv.data(), nbf, nbf};

        // fock = Hcore + Juv + Kuv; fock = (fock + fock.transpose()) * 0.5;
        fock = Hcore + Juv - Kuv; fock = (fock + fock.transpose()) * 0.5;
        focks.emplace_back(fock);
        diis_res.emplace_back(sij_inv_half * (fock * Puv * sij_map - sij_map * Puv * fock) * sij_inv_half);
        auto adaptor_fock = xt::adapt(focks.back().data(), nbf*nbf, xt::no_ownership(), xt::xtensor<double, 2>::shape_type{static_cast<size_t>(nbf), static_cast<size_t>(nbf)});
        auto adaptor_res  = xt::adapt(diis_res.back().data(), nbf*nbf, xt::no_ownership(), xt::xtensor<double, 2>::shape_type{static_cast<size_t>(nbf), static_cast<size_t>(nbf)});

        double kin_e = (Puv * tij_map.transpose()).trace();
        double ext_e = (Puv * vij_map.transpose()).trace();
        double hartree_e  = (Puv * Juv.transpose()).trace() * 0.5;
        // double exchange_e = (Puv * Kuv.transpose()).trace() * 0.5;
        double exchange_e = (Puv * Kuv.transpose()).trace() * -0.5;
        double etot = 0.5 * (Puv.array() * (Hcore + fock).array()).sum() + e_nuc; // (A * B).sum() <=> (A @ B.T).trace()

        double e_diff = etot - etot_old;
        double d_diff = std::sqrt((Puv.array() - Puv_old.array()).square().sum());

        std::cout << "One-Electron energy = " << std::setw(25) << std::setprecision(15) << std::fixed << std::right << kin_e + ext_e          << std::endl;
        std::cout << "Two-Electron energy = " << std::setw(25) << std::setprecision(15) << std::fixed << std::right << hartree_e + exchange_e << std::endl;
        std::cout << "Total energy        = " << std::setw(25) << std::setprecision(15) << std::fixed << std::right << etot                   << std::endl;
        // std::cout << std::setprecision(15) << kin_e + ext_e + hartree_e + exchange_e + e_nuc << std::endl;
        std::cout << "Energy difference   = " << std::setw(25) << std::setprecision(15) << std::fixed << std::right << e_diff                 << std::endl;
        std::cout << "Density difference  = " << std::setw(25) << std::setprecision(15) << std::fixed << std::right << d_diff                 << std::endl;

        // DIIS
        if (counter > 2) diis(fock, focks, diis_res);


#ifdef RELEASE
        if (std::abs(etot - etot_old) <= e_convergence && std::abs(d_diff) <= d_convergence) {
            std::cout << std::endl;
            std::cout << "Doubly occupied:\n";
            int iorb = 1;
            for (auto val : vals.head(nocc).transpose()) {
                std::cout << std::setw(15) << std::setprecision(8) << std::fixed << std::right << val << "    ";
                if (iorb % 4 == 0) std::cout << std::endl;
                iorb += 1;
            }

            iorb = 1;
            std::cout << "\nVirtual:\n";
            for (auto val : vals.tail(nbf-nocc).transpose()) {
                std::cout << std::setw(15) << std::setprecision(8) << std::fixed << std::right << val << "    ";
                if (iorb % 4 == 0) std::cout << std::endl;
                iorb += 1;
            }
            std::cout << "\n>! NORMAL TERMINALTION\n";
            break;
        }

        counter += 1;
        etot_old = etot;
        Puv_old = Puv;

        if (counter > maxiter) {std::cout << "MAXITER exceeded\n"; break;}
    }
#endif


    std::cout << ">! safe here\n";
}


