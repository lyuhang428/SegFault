#include "../include/mints.hpp"

double sf::sfact2(int n)
{
    if (n <= 0) return 1.;
    return n * sfact2(n - 2);
}

std::vector<std::array<int, 3>> sf::cart_ordering(int ltot)
{

    assert(ltot >= 0);
    int nbf_cart = (ltot + 1) * (ltot + 2) / 2;
    // int nbf_pure = 2 * l + 1;

    std::vector<std::array<int, 3>> table_cart;
    table_cart.reserve(nbf_cart);

    for (auto lx=ltot; lx >= 0; --lx) {
        for (auto ly=ltot - lx; ly >= 0; --ly) {
            auto lz = ltot - lx - ly;
            // std::string x(lx, 'x');
            // std::string y(ly, 'y');
            // std::string z(lz, 'z');
            // std::cout << "(" << lx << "," << ly << "," << lz << ")\n";
            // std::cout << x << y << z << "\n";
            table_cart.push_back({lx, ly, lz});
        }
    }

    return table_cart;
}


std::vector<std::vector<libint2::Shell>> sf::read_g94_basis_library(std::string file_dot_g94,
                                                                           bool force_cartesian_d,
                                                                           bool throw_if_missing,
                                                                           std::string locale_name)
{
    std::locale locale(locale_name.c_str());
    std::vector<std::vector<libint2::Shell>> ref_shells(119); // 118 = number of chemical elements + 1 because we insert at(Z) // TO BE RETURNED
    std::ifstream is(file_dot_g94);
    is.imbue(locale);

    if (is.good()) {
        std::string line, rest;

        auto LIBINT2_LINE_TO_STRINGSTREAM = [&](const std::string& line)
        {
            std::istringstream iss{line};
            iss.imbue(locale);
            return iss;
        };

        size_t Z;
        bool nextbasis = true;
        bool nextshell = false;
        bool first_element = true;

        // read lines till end
        do {
            // skipping empties and starting with '!' (the comment delimiter)
            if (line.empty() || line[0] == '!') continue;

            if (line == "****") {
                // old (EMSL) basis set exchange g94 format marks the beginning of data by ****
                // new (MolSSI) basis set exchange g94 format does not start with ****
                // so if found **** and still waiting for the first element, this is new g94 format, skip to next line
                if (first_element) continue;
                nextbasis = true;
                nextshell = false;
                continue;
            }

            if (nextbasis) {
                nextbasis = false;
                first_element = false;
                auto iss = LIBINT2_LINE_TO_STRINGSTREAM(line);
                std::string elemsymbol;
                iss >> elemsymbol >> rest;

                bool found = false;
                for (const auto &e: libint2::chemistry::get_element_info()) {
                    if (strcaseequal(e.symbol, elemsymbol)) {
                        Z = e.Z;
                        found = true;
                        break;
                    }
                }

                if (not found) {
                    std::ostringstream oss;
                    oss << "in file " << file_dot_g94 << " found G94 basis set for element symbol \"" << elemsymbol << "\", not found in Periodic Table.";
                    throw std::logic_error(oss.str());
                }

                nextshell = true;
                continue;
            }

            if (nextshell) {
                auto iss = LIBINT2_LINE_TO_STRINGSTREAM(line);
                std::string amlabel;
                std::size_t nprim;
                iss >> amlabel >> nprim >> rest;
                if (amlabel != "SP" && amlabel != "sp") {
                    assert(amlabel.size() == 1);
                    // Gaussian labels L=6,7,8,9... AOs as I,J,K,L... instead of I,K,L,M..., see https://github.com/MolSSI-BSE/basis_set_exchange/issues/292
                    const auto amlabel0_is_J = amlabel[0] == 'j' || amlabel[0] == 'J';
                    auto l = (amlabel0_is_J) ? 7 : libint2::Shell::am_symbol_to_l(amlabel[0]);
                    if (l >= 7 && !amlabel0_is_J) l++;  // Gaussian's 'K' means 'L', 'L' means 'M', etc.
                    std::vector<double> exps;
                    std::vector<double> coeffs;
                    for (decltype(nprim) p = 0; p != nprim; ++p) {
                        while (std::getline(is, line) && (line.empty() || line[0] == '!')) {}
                        auto iss = LIBINT2_LINE_TO_STRINGSTREAM(line);
                        double e, c;
                        iss >> e >> c;
                        exps.emplace_back(e);
                        coeffs.emplace_back(c);
                    }
                    auto pure = force_cartesian_d ? (l > 2) : (l > 1);
                    ref_shells.at(Z).push_back(libint2::Shell{std::move(exps), {{l, pure, std::move(coeffs)}}, {{0, 0, 0}}, false}); // false!
                }
                else { // split the SP shells
                    std::vector<double> exps;
                    std::vector<double> coeffs_s, coeffs_p;
                    for (decltype(nprim) p = 0; p != nprim; ++p) {
                        while (std::getline(is, line) && (line.empty() || line[0] == '!')) {}
                        auto iss = LIBINT2_LINE_TO_STRINGSTREAM(line);
                        double e, c1, c2;
                        iss >> e >> c1 >> c2;
                        exps.emplace_back(e);
                        coeffs_s.emplace_back(c1);
                        coeffs_p.emplace_back(c2);
                    }
                    // if false, DO NOT embed normalization into coefficient
                    // if false, return coefficient as it is in .g94 file
                    // embedding will be done OUTSIDE this fn when manually build libint2::Shell
                    ref_shells.at(Z).push_back(libint2::Shell{exps, {{0, false, coeffs_s}}, {{0, 0, 0}}, false});
                    ref_shells.at(Z).push_back(libint2::Shell{std::move(exps), {{1, false, std::move(coeffs_p)}}, {{0, 0, 0}}, false});
                }
            }
        } while (std::getline(is, line));

    }
    else {  // !is.good()
        if (throw_if_missing) {
            std::ostringstream oss;
            oss << "BasisSet::read_g94_basis_library(): could not open \"" << file_dot_g94 << "\"" << std::endl;
            throw std::ios_base::failure(oss.str());
        }
    }

    return ref_shells;
}

std::vector<libint2::Atom> sf::make_atoms(const std::string& xyzfile)
{
    // .xyz angstrom automatically converted to Bohr
    std::ifstream fp_mol(xyzfile, std::ios::in);
    assert(fp_mol.is_open());
    auto atoms = libint2::read_dotxyz(fp_mol);
    return atoms;
}

std::vector<libint2::Shell> sf::make_shells(const std::vector<libint2::Atom>& atoms, bool pure, const std::string& file_dot_g94)
{
    // 提取对应的指数 系数 角动量
    // 通过该方式得到的所有 Shell 坐标都在原点
    // 某元素的所有 Shell
    // 该函数内部设置不将归一化常数合并如缩并系数中，基组文件什么样返回的就是什么样
    // 对于 SP 类型的数据，read_g94_basis_library 已将其拆分
    // 通过该函数返回的 Shell::contr[0].coeff 已经混合了归一化常数，不能用于基函数在网格上的计算
    std::vector<std::vector<libint2::Shell>> basis = read_g94_basis_library(file_dot_g94); // 读取整个基组文件，若某一元素的基组信息不包含在内，其长度为0

    std::vector<libint2::Shell> shells; // TO BE RETURNED

    for (const auto& atom : atoms) {
        assert(basis[atom.atomic_number].size() > 0); // check if element exists in the given basis file 不存在的话长度为0
        // loop all shells in an atom, extract exponents, coefficients, l, xyz, and enable embedding normalization into coefficients
        for (const auto& shell : basis[atom.atomic_number]) {
            assert(shell.contr.size() == 1); // SP 类型的数据应该已被拆分
            shells.emplace_back(libint2::Shell{shell.alpha, {{shell.contr[0].l, pure, shell.contr[0].coeff}}, std::array<double, 3>{atom.x, atom.y, atom.z}, true}); // if true, embed normalization into coefficient
        }
    }
    return shells;
}

std::vector<libint2::Shell> sf::_make_shells_cart_noembed(const std::vector<libint2::Atom>& atoms, const std::string& file_dot_g94)
{
    std::vector<std::vector<libint2::Shell>> basis = read_g94_basis_library(file_dot_g94); // 读取整个基组文件，若某一元素的基组信息不包含在内，其长度为0

    std::vector<libint2::Shell> shells; // TO BE RETURNED

    for (const auto& atom : atoms) {
        assert(basis[atom.atomic_number].size() > 0); // check if element exists in the given basis file 不存在的话长度为0
        // loop all shells in an atom, extract exponents, coefficients, l, xyz, and enable embedding normalization into coefficients
        for (const auto& shell : basis[atom.atomic_number]) {
            assert(shell.contr.size() == 1); // SP 类型的数据应该已被拆分
            shells.emplace_back(libint2::Shell{shell.alpha, {{shell.contr[0].l, false, shell.contr[0].coeff}}, std::array<double, 3>{atom.x, atom.y, atom.z}, false}); // if true, embed normalization into coefficient
        }
    }
    return shells;
}

std::vector<size_t> sf::get_shell2bf(const std::vector<libint2::Shell>& shells)
{
    std::vector<size_t> shell2bf; // to be returned
    shell2bf.reserve(shells.size());
    int id0 = 0;
    for (const auto& shell : shells) {
        shell2bf.emplace_back(id0);
        id0 += shell.size();
    }
    return shell2bf;
}

size_t sf::get_max_nprim(const std::vector<libint2::Shell>& shells)
{
    size_t max_nprim = 0;
    for (const auto& shell : shells) max_nprim > shell.nprim() ?  : max_nprim = shell.nprim();
    assert(max_nprim > 0); // at least 1
    return max_nprim;
}

int sf::get_lmax(const std::vector<libint2::Shell>& shells)
{
    int lmax = -1;
    for (const auto& shell : shells) {
        // 通常来说 contr 长度为1,除非是 SP 类型
        for (const auto& contr : shell.contr) lmax > contr.l ? : lmax = contr.l;
    }
    assert(lmax >= 0); // cannot be negative
    return lmax;
}

xt::xtensor<double, 2> sf::get_olp(const std::vector<libint2::Shell>& shells)
{
    int nshell = shells.size(); // number of shell ; i.e. (SSSPPD) = 6
    int nbf = 0;
    for (const auto& shell : shells) nbf += shell.size(); // number of basis function ; i.e. (SSSPPD)_cart = 15 or (SSSPPD)_pure = 14

    xt::xtensor<double, 2> sij = xt::zeros<double>({nbf, nbf}); // to be returned

    const std::vector<size_t> shell2bf = get_shell2bf(shells);  // bf index
    const size_t             max_nprim = get_max_nprim(shells); // 最大缩并 ; i.e. S 8 1.00 / S 1 1.00 then max_nprim = 8
    const int                     lmax = get_lmax(shells);      // 最高角动量

    auto engine = libint2::Engine{libint2::Operator::overlap, max_nprim, lmax, 0}; // 0 means no deriv
    engine.set(libint2::CartesianShellNormalization::uniform);

    // for olp len(buf_vec)=0 and buf_vec[0] is double*, loop over buf_vec[0] get integral values
    // if for dipole integral, len(buf_vec)=4, [0] is olp, [1] is <|x|>, [2] is <|y|>, etc.
    // if for quadruple integral, len=10, [0] is olp, [1-3] is dipole, [4-] is xx, xy, xz, yy, yz, zz
    // target_ptr_vec == std::vector<const value_type*, detail::ext_stack_allocator<const value_type*, max_ntargets>>;
    // where value_type = LIBINT2_REALTYPE; numeric.h
    // where #define LIBINT2_REALTYPE double; config.h
    // so buf_vec is std::vector<const double*>, and its length for olp, kin, ext, eri is 1
    // double* points and stores all integral values computed
    const auto& buf_vec = engine.results();

    for (auto s1=0; s1 < nshell; ++s1) {
        for (auto s2=s1; s2 < nshell; ++s2) {
            engine.compute(shells[s1], shells[s2]); // 按照 shell 进行块计算 block
            const double* ints = buf_vec[0]; // 1D double array, need to reshape it in row-major (C)
            if (ints == nullptr) continue;
            const size_t ibf1 = shell2bf[s1];      // bf id
            const size_t ibf2 = shell2bf[s2];      // bf id
            const size_t nbf1 = shells[s1].size(); // number of bf in this shell
            const size_t nbf2 = shells[s2].size(); // number of bf in this shell
            const auto adaptor = xt::adapt(ints, nbf1*nbf2, xt::no_ownership(), xt::xtensor<double, 2>::shape_type{nbf1, nbf2});
            xt::view(sij, xt::range(ibf1, ibf1+nbf1), xt::range(ibf2, ibf2+nbf2)) = adaptor; // copy double* to xtensor
            if (s2 != s1) xt::view(sij, xt::range(ibf2, ibf2+nbf2), xt::range(ibf1, ibf1+nbf1)) = xt::transpose(adaptor);
        }
    }

    return sij;
}

xt::xtensor<double, 2> sf::get_kin(const std::vector<libint2::Shell>& shells)
{
    int nshell = shells.size();
    int nbf = 0;
    for (const auto& shell : shells) nbf += shell.size();

    xt::xtensor<double, 2> tij = xt::zeros<double>({nbf, nbf}); // to be returned

    const std::vector<size_t> shell2bf = get_shell2bf(shells);
    const size_t             max_nprim = get_max_nprim(shells);
    const int                     lmax = get_lmax(shells);

    auto engine = libint2::Engine{libint2::Operator::kinetic, max_nprim, lmax, 0};
    engine.set(libint2::CartesianShellNormalization::uniform);
    const auto& buf_vec = engine.results();

    for (auto s1=0; s1 < nshell; ++s1) {
        for (auto s2=s1; s2 < nshell; ++s2) {
            engine.compute(shells[s1], shells[s2]);
            const double* ints = buf_vec[0];
            if (ints == nullptr) continue;
            size_t ibf1 = shell2bf[s1];
            size_t ibf2 = shell2bf[s2];
            size_t nbf1 = shells[s1].size();
            size_t nbf2 = shells[s2].size();
            const auto adaptor = xt::adapt(ints, nbf1*nbf2, xt::no_ownership(), xt::xtensor<double, 2>::shape_type{nbf1, nbf2});
            xt::view(tij, xt::range(ibf1, ibf1+nbf1), xt::range(ibf2, ibf2+nbf2)) = adaptor;
            if (s2 != s1) xt::view(tij, xt::range(ibf2, ibf2+nbf2), xt::range(ibf1, ibf1+nbf1)) = xt::transpose(adaptor);
        }
    }

    return tij;
}

xt::xtensor<double, 2> sf::get_ext(const std::vector<libint2::Shell>& shells, const std::vector<libint2::Atom>& atoms)
{
    int nshell = shells.size();
    int nbf = 0;
    for (const auto&shell : shells) nbf += shell.size();

    xt::xtensor<double, 2> vij = xt::zeros<double>({nbf, nbf});

    const std::vector<size_t> shell2bf = get_shell2bf(shells);
    const size_t             max_nprim = get_max_nprim(shells);
    const int                     lmax = get_lmax(shells);

    auto engine = libint2::Engine{libint2::Operator::nuclear, max_nprim, lmax, 0};
    engine.set(libint2::CartesianShellNormalization::uniform);
    engine.set_params(libint2::make_point_charges(atoms));
    const auto& buf_vec = engine.results();

    for (auto s1=0; s1 < nshell; ++s1) {
        for (auto s2=s1; s2 < nshell; ++s2) {
            engine.compute(shells[s1], shells[s2]);
            const double* ints = buf_vec[0];
            if (ints == nullptr) continue;
            size_t ibf1 = shell2bf[s1];
            size_t ibf2 = shell2bf[s2];
            size_t nbf1 = shells[s1].size();
            size_t nbf2 = shells[s2].size();
            const auto adaptor = xt::adapt(ints, nbf1*nbf2, xt::no_ownership(), xt::xtensor<double, 2>::shape_type{nbf1, nbf2});
            xt::view(vij, xt::range(ibf1, ibf1+nbf1), xt::range(ibf2, ibf2+nbf2)) = adaptor;
            if (s2 != s1) xt::view(vij, xt::range(ibf2, ibf2+nbf2), xt::range(ibf1, ibf1+nbf1)) = xt::transpose(adaptor);
        }
    }

    return vij;
}

std::pair<xt::xtensor<double, 2>, xt::xtensor<double, 2>> sf::get_JK(const std::vector<libint2::Shell>& shells, const xxd& D)
{
    const int nshell = shells.size();
    int nbf = 0;
    for (const auto& shell : shells) nbf += shell.size();
    assert(D.rows() == nbf && D.cols() == nbf);
    const int             max_nprim = get_max_nprim(shells);
    const int                  lmax = get_lmax(shells);
    const std::vector<size_t> shell2bf = get_shell2bf(shells);

    xt::xtensor<double, 2> Juv = xt::zeros<double>({nbf, nbf}); // to be returned
    xt::xtensor<double, 2> Kuv = xt::zeros<double>({nbf, nbf}); // to be returned

    auto engine = libint2::Engine(libint2::Operator::coulomb, max_nprim, lmax);
    engine.set(libint2::CartesianShellNormalization::uniform);
    const auto& buf = engine.results();


    // loop over permutationally-unique set of shells
    for (auto s1 = 0; s1 != nshell; ++s1) {
        const auto bf1_first = shell2bf[s1];  // first basis function in this shell
        const auto n1 = shells[s1].size();    // number of basis functions in this shell
        for (auto s2 = 0; s2 <= s1; ++s2) {
            const auto bf2_first = shell2bf[s2];
            const auto n2 = shells[s2].size();
            for (auto s3 = 0; s3 <= s1; ++s3) {
                const auto bf3_first = shell2bf[s3];
                const auto n3 = shells[s3].size();
                const auto s4_max = (s1 == s3) ? s2 : s3;
                for (auto s4 = 0; s4 <= s4_max; ++s4) {
                    const auto bf4_first = shell2bf[s4];
                    const auto n4 = shells[s4].size();
                    // compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
                    const auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
                    const auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
                    const auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
                    const auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

                    // engine.compute is not thread-safe
                    // to run it in parallel, each thread has to make a new engine and do compute
                    engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
                    const auto* buf_1234 = buf[0];
                    if (buf_1234 == nullptr) continue;  // if all integrals screened out, skip to next quartet

                    for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
                        const auto bf1 = f1 + bf1_first;
                        for (auto f2 = 0; f2 != n2; ++f2) {
                            const auto bf2 = f2 + bf2_first;
                            for (auto f3 = 0; f3 != n3; ++f3) {
                                const auto bf3 = f3 + bf3_first;
                                for (auto f4 = 0; f4 != n4; ++f4, ++f1234) {
                                    const auto bf4               = f4 + bf4_first;
                                    const auto value             = buf_1234[f1234]; // 积分结果
                                    const auto value_scal_by_deg = value * s1234_deg;
                                    // G(bf1, bf2) += 0.5   * D(bf3, bf4) * value_scal_by_deg;
                                    // G(bf3, bf4) += 0.5   * D(bf1, bf2) * value_scal_by_deg;
                                    // G(bf1, bf3) -= 0.125 * D(bf2, bf4) * value_scal_by_deg;
                                    // G(bf2, bf4) -= 0.125 * D(bf1, bf3) * value_scal_by_deg;
                                    // G(bf1, bf4) -= 0.125 * D(bf2, bf3) * value_scal_by_deg;
                                    // G(bf2, bf3) -= 0.125 * D(bf1, bf4) * value_scal_by_deg;
                                    Juv(bf1, bf2) += 0.5   * D(bf3, bf4) * value_scal_by_deg;
                                    Juv(bf3, bf4) += 0.5   * D(bf1, bf2) * value_scal_by_deg;
                                    // Kuv(bf1, bf3) -= 0.125 * D(bf2, bf4) * value_scal_by_deg;
                                    // Kuv(bf2, bf4) -= 0.125 * D(bf1, bf3) * value_scal_by_deg;
                                    // Kuv(bf1, bf4) -= 0.125 * D(bf2, bf3) * value_scal_by_deg;
                                    // Kuv(bf2, bf3) -= 0.125 * D(bf1, bf4) * value_scal_by_deg;
                                    Kuv(bf1, bf3) += 0.125 * D(bf2, bf4) * value_scal_by_deg;
                                    Kuv(bf2, bf4) += 0.125 * D(bf1, bf3) * value_scal_by_deg;
                                    Kuv(bf1, bf4) += 0.125 * D(bf2, bf3) * value_scal_by_deg;
                                    Kuv(bf2, bf3) += 0.125 * D(bf1, bf4) * value_scal_by_deg;
                                }
                            }
                        }
                    } // end value assignment
                }
            }
        }
    }
    Juv = (Juv + xt::transpose(Juv)) * 0.5;
    Kuv = (Kuv + xt::transpose(Kuv)) * 0.5;
    return {Juv, Kuv};
}



sf::Molecule::Molecule(const std::string& xyzfile, const std::string& name) : xyzfile(xyzfile), name(name) {
    this->atoms = make_atoms(xyzfile);
    this->shells_pure = make_shells(this->atoms, true, this->name);
    this->shells_cart = make_shells(this->atoms, false, this->name);
    this->shells_cart_noembed = _make_shells_cart_noembed(this->atoms, this->name);
    assert(this->shells_cart.size() == this->shells_pure.size()); // 无论球谐基函数还是笛卡尔基函数都应该有相同的 nshell

    this->nshell = shells_pure.size();
    this->natom = this->atoms.size();

    this->xyz.resize({static_cast<size_t>(this->natom), 3});
    this->ne = 0;

    for (auto iatom=0; iatom < this->natom; ++iatom) {
        this->numbers.emplace_back(this->atoms[iatom].atomic_number);
        this->ne += this->atoms[iatom].atomic_number;
        xt::row(this->xyz, iatom) = xt::xtensor_fixed<double, xt::xshape<3>>{this->atoms[iatom].x, this->atoms[iatom].y, this->atoms[iatom].z};
        this->symbols.emplace_back(::Number2Atom.at(this->atoms[iatom].atomic_number));

        if (this->atoms[iatom].atomic_number == 1) this->rms.emplace_back(::BraggSlaterRadii.at("H") * ::ANG2BOHR);
        else this->rms.emplace_back(::BraggSlaterRadii.at(this->symbols[iatom]) * 0.5 * ::ANG2BOHR);
    }

    this->nocc = this->ne / 2;

    this->nbf_cart = 0;
    this->nbf_pure = 0;
    for (const auto& shell : this->shells_cart) {
        this->nbf_cart += shell.size();
        this->nbf_in_shells_cart.emplace_back(shell.size()); // for building cart2pure tf
    }
    for (const auto& shell : this->shells_pure) {
        this->nbf_pure += shell.size();
        this->nbf_in_shells_pure.emplace_back(shell.size()); // for building cart2pure tf
    }

    this->e_nuc = get_e_nuc();
}

xxd sf::Molecule::make_tf() const
{
    const int nrow = this->nbf_pure;
    const int ncol = this->nbf_cart;
    xxd TF = xxd::Zero(nrow, ncol); // to be returned 块对角矩阵，高度稀疏

    int row_offset = 0;
    int col_offset = 0;

    for (auto ishell=0; ishell < this->nshell; ++ishell) {
        if (this->nbf_in_shells_pure[ishell] == 1 && this->nbf_in_shells_cart[ishell] == 1) { // s-orbital
            TF(row_offset, col_offset) = 1.;
            row_offset++;
            col_offset++;
        }
        else if (this->nbf_in_shells_pure[ishell] == 3 && this->nbf_in_shells_cart[ishell] == 3) { // p
            TF.block(row_offset, col_offset, 3, 3) = tf1;
            row_offset += 3;
            col_offset += 3;
        }
        else if (this->nbf_in_shells_pure[ishell] == 5 && this->nbf_in_shells_cart[ishell] == 6) { // d
            TF.block(row_offset, col_offset, 5, 6) = tf2;
            row_offset += 5;
            col_offset += 6;
        }
        else if (this->nbf_in_shells_pure[ishell] == 7 && this->nbf_in_shells_cart[ishell] == 10) { // f
            TF.block(row_offset, col_offset, 7, 10) = tf3;
            row_offset += 7;
            col_offset += 10;
        }
        else {std::cerr << "Beyond f-orbital not available\n"; exit(-1);};
    }

    return TF;
}

xt::xtensor<double, 2> sf::Molecule::make_tf_xt() const
{
    const int nrow = this->nbf_pure;
    const int ncol = this->nbf_cart;
    xt::xtensor<double, 2> TF = xt::zeros<double>({nrow, ncol}); // to be returned 块对角矩阵，高度稀疏

    int row_offset = 0;
    int col_offset = 0;

    for (auto ishell=0; ishell < this->nshell; ++ishell) {
        if (this->nbf_in_shells_pure[ishell] == 1 && this->nbf_in_shells_cart[ishell] == 1) { // s-orbital
            TF(row_offset, col_offset) = 1.;
            row_offset++;
            col_offset++;
        }
        else if (this->nbf_in_shells_pure[ishell] == 3 && this->nbf_in_shells_cart[ishell] == 3) { // p
            xt::view(TF, xt::range(row_offset, row_offset+3), xt::range(col_offset, col_offset+3)) = xt::adapt(tf1.data(), tf1.size(), xt::no_ownership(), xt::xtensor<double, 2>::shape_type{3,3});
            row_offset += 3;
            col_offset += 3;
        }
        else if (this->nbf_in_shells_pure[ishell] == 5 && this->nbf_in_shells_cart[ishell] == 6) { // d
            xt::view(TF, xt::range(row_offset, row_offset+5), xt::range(col_offset, col_offset+6)) = xt::adapt(tf2.data(), tf2.size(), xt::no_ownership(), xt::xtensor<double, 2>::shape_type{5,6});
            row_offset += 5;
            col_offset += 6;
        }
        else if (this->nbf_in_shells_pure[ishell] == 7 && this->nbf_in_shells_cart[ishell] == 10) { // f
            xt::view(TF, xt::range(row_offset, row_offset+7), xt::range(col_offset, col_offset+10)) = xt::adapt(tf3.data(), tf3.size(), xt::no_ownership(), xt::xtensor<double, 2>::shape_type{7,10});
            row_offset += 7;
            col_offset += 10;
        }
        else {std::cerr << "Beyond f-orbital not available\n"; exit(-1);};
    }

    return TF;
}

double sf::Molecule::get_e_nuc() const
{
    double _e_nuc = 0;
    for (auto iatom=0; iatom < this->atoms.size()-1; ++iatom) {
        for (auto jatom=iatom+1; jatom < this->atoms.size(); ++jatom) {
            xt::xtensor_fixed<double, xt::xshape<3>> ri{this->atoms[iatom].x, this->atoms[iatom].y, this->atoms[iatom].z};
            xt::xtensor_fixed<double, xt::xshape<3>> rj{this->atoms[jatom].x, this->atoms[jatom].y, this->atoms[jatom].z};
            _e_nuc += this->atoms[iatom].atomic_number * this->atoms[jatom].atomic_number / xt::norm_l2(ri - rj)();
        }
    }
    return _e_nuc;
}



sf::BFs::BFs(const sf::Molecule& mol)
{
    std::vector<BF> _bfs; // to be moved
    for (const auto& shell : mol.shells_cart_noembed) {
        // shells_cart_noembed 仅用来记录基组文件信息，赋值完成后不用于任何计算
        int l = shell.contr[0].l; // total l, lx, ly, lz to be given by cart_ordering
        std::vector<std::array<int, 3>> orderings = cart_ordering(l); // i.e. s x y z xx xy xz yy yz zz xxx xxy xxz xyy ...
        xt::xtensor<double, 1> exponents    = xt::adapt(shell.alpha.data(), shell.alpha.size(), xt::no_ownership(), xt::xtensor<double, 1>::shape_type{shell.alpha.size()});
        xt::xtensor<double, 1> coefficients = xt::adapt(shell.contr[0].coeff.data(), shell.contr[0].coeff.size(), xt::no_ownership(), xt::xtensor<double, 1>::shape_type{shell.contr[0].coeff.size()});
        xt::xtensor_fixed<double, xt::xshape<3>> center{shell.O[0], shell.O[1], shell.O[2]}; // in Bohr automatically
        auto _a = xt::eval(xt::pow(2. * exponents / PI, 0.75));
        auto _b = xt::eval(xt::pow(4. * exponents, l * 0.5));
        for (const auto& ordering : orderings) {
            double _c = std::sqrt(sfact2(2 * ordering[0] - 1) * sfact2(2 * ordering[1] - 1) * sfact2(2 * ordering[2] - 1));
            xt::xtensor<double, 1> nrfcs = _a * _b / _c; // libint2::Engine.set(Engine::CartesianShellNormalization::uniform) to be consistent
            _bfs.emplace_back(BF{exponents, coefficients, nrfcs, center, xt::xtensor_fixed<int, xt::xshape<3>>{ordering[0], ordering[1], ordering[2]}});
        }
    }
    this->bfs = std::move(_bfs);
    this->nbf = this->bfs.size();
}

xt::xtensor<double, 1> sf::BFs::ao_val(double x, double y, double z) const
{
    xt::xtensor_fixed<double, xt::xshape<3>> xyz = {x, y, z};
    xt::xtensor<double, 1> res_cart = xt::zeros<double>({static_cast<size_t>(this->nbf)}); // to be returned {nbf, }

    size_t ibf = 0;
    for (const auto& bf : this->bfs) {
        double _angular = xt::prod(xt::pow((xyz - bf.center), bf.lxyz))();
        xt::xtensor<double, 1> _radial = xt::exp(-1. * bf.exponents * xt::sum(xt::square(xyz - bf.center)));
        res_cart[ibf] = xt::sum(bf.coefficients * bf.normfactors * _angular * _radial)();
        ibf++;
    }

    return res_cart;
}

xt::xtensor<double, 2> sf::BFs::ao_val(const xt::xtensor<double, 1>& x, const xt::xtensor<double, 1>& y, const xt::xtensor<double, 1>& z) const
{
    assert(x.size() == y.size() && x.size() == z.size());
    size_t npoints = x.size();

    xt::xtensor<double, 2> res_cart = xt::zeros<double>({static_cast<size_t>(this->nbf), npoints}); // to be returned {nbf, npoints}

    size_t ibf = 0;
    for (const auto& bf : this->bfs) {
        xt::xtensor<double, 1> _angular = xt::pow(x - bf.center[0], bf.lxyz[0]) * xt::pow(y - bf.center[1], bf.lxyz[1]) * xt::pow(z - bf.center[2], bf.lxyz[2]);
        for (auto iprim=0; iprim < bf.nprim; ++iprim) {
            xt::xtensor<double, 1> _radial = xt::exp(-bf.exponents[iprim] * (xt::square(x - bf.center[0]) + xt::square(y - bf.center[1]) + xt::square(z - bf.center[2])));
            xt::row(res_cart, ibf) += bf.coefficients[iprim] * bf.normfactors[iprim] * _angular * _radial;
        }
        ibf++;
    }

    return res_cart;
}

