#include "../include/beckegrid.hpp"


double beckegrid::sk(const int k, const double nuij)
{
    assert(1 <= k && k <= 6);
    switch (k)
    {
    case 1:
        return beckegrid::s1(nuij);
    case 2:
        return beckegrid::s2(nuij);
    case 3:
        return beckegrid::s3(nuij);
    case 4:
        return beckegrid::s4(nuij);
    case 5:
        return beckegrid::s5(nuij);
    case 6:
        return beckegrid::s6(nuij);
    default:
        exit(-1);
    }
}

void beckegrid::sk(const int k, const xt::xtensor<double, 1>& nuij, xt::xtensor<double, 1>& sij)
{
    assert(1 <= k && k <= 6);
    switch (k)
    {
    case 1:
        std::transform(nuij.begin(), nuij.end(), sij.begin(), beckegrid::s1);
        break;
    case 2:
        std::transform(nuij.begin(), nuij.end(), sij.begin(), beckegrid::s2);
        break;
    case 3:
        std::transform(nuij.begin(), nuij.end(), sij.begin(), beckegrid::s3);
        break;
    case 4:
        std::transform(nuij.begin(), nuij.end(), sij.begin(), beckegrid::s4);
        break;
    case 5:
        std::transform(nuij.begin(), nuij.end(), sij.begin(), beckegrid::s5);
        break;
    case 6:
        std::transform(nuij.begin(), nuij.end(), sij.begin(), beckegrid::s6);
        break;
    default:
        exit(-1);
    }
}

xt::xtensor<double, 2> beckegrid::gaussCheby2(const size_t norder, const double rm)
{
    xt::xtensor<double, 2> zxrw{{4, norder}, 0.}; // to be returned
    xt::row(zxrw, 0) = xt::arange(1., static_cast<double>(norder+1)); // z = 1,2,3,...,norder
    xt::row(zxrw, 1) = xt::cos(xt::row(zxrw, 0) / (norder + 1) * PI);      // x
    xt::row(zxrw, 2) = rm * (1. + xt::row(zxrw, 1)) / (1. - xt::row(zxrw, 1));  // r
    xt::row(zxrw, 3) = PI / (norder + 1) * xt::square(xt::sin(xt::row(zxrw, 0) / (norder + 1) * PI)) / xt::sqrt(1. - xt::square(xt::row(zxrw, 1))) * 2. * rm / xt::square(1. - xt::row(zxrw, 1)); // w

    return zxrw;
}

beckegrid::BeckeFuzzyCell::BeckeFuzzyCell(const std::vector<std::string>& symbols,
                                          const std::vector<double>&          rms,
                                          const std::vector<int>&         numbers,
                                          const xt::xtensor<double, 2>&       xyz,
                                          const size_t                       nrad,
                                          const size_t                       nleb,
                                          const size_t                          k,
                                          const bool                       biased)
: symbols(symbols),
          rms(rms),
  numbers(numbers),
          xyz(xyz),
        nrad(nrad),
        nleb(nleb),
              k(k),
     biased(biased)
{
    assert((symbols.size() == rms.size()) && symbols.size() == numbers.size() && symbols.size() == xyz.shape(0));
    assert(xyz.shape(1) == 3);
    assert(k>=1 && k<=6);
    this->natom = symbols.size();
}

xt::xtensor<double, 1> beckegrid::BeckeFuzzyCell::get_weight_s(const double x, const double y, const double z)
{
    xt::xtensor_fixed<double, xt::xshape<3>> r{x, y, z};
    xt::xtensor<double, 1> pres = xt::zeros<double>({natom}); // to be returned

    for (auto iatom=0; iatom < this->natom; ++iatom)
    {
        double s = 1.;
        for (auto jatom=0; jatom < this->natom; ++jatom)
        {
            if (jatom == iatom) continue;
            double ri   = xt::norm_l2(r - xt::row(xyz, iatom))();
            double rj   = xt::norm_l2(r - xt::row(xyz, jatom))();
            double rij  = xt::norm_l2(xt::row(xyz, iatom) - xt::row(xyz, jatom))();
            double muij = (ri - rj) / rij;
            double nuij = muij;
            double sij  = 0.;

            if (this->biased)
            {
                double chi = BraggSlaterRadii.at(this->symbols[iatom]) / BraggSlaterRadii.at(this->symbols[jatom]);
                double uij = (chi - 1.) / (chi + 1.);
                double aij = uij / (uij * uij - 1);
                if (aij > 0.5) aij = 0.5;
                else if (aij < -0.5) aij = -0.5;
                nuij = muij + aij * (1. - muij * muij);
            }

            sij = beckegrid::sk(this->k, nuij);
            s *= sij;
        }
        pres[iatom] = s;
    }
    pres /= xt::sum(pres, {0});
    return pres;
}

xt::xtensor<double, 2> beckegrid::BeckeFuzzyCell::get_weight_p(const xt::xtensor<double, 2>& grid)
{
    assert(grid.shape(1) == 3); // (N, 3)
    xt::xtensor<double, 2> dist = xt::norm_l2(this->xyz - xt::view(grid, xt::all(), xt::newaxis(), xt::all()), {2}); // (ngrid, natom)
    xt::xtensor<double, 2> pres = xt::zeros<double>({this->natom, grid.shape(0)}); // to be returned (natom, ngrid)

    for (auto iatom=0; iatom < this->natom; ++iatom)
    {
        xt::xtensor<double, 1> s = xt::ones<double>({grid.shape(0)});
        for (auto jatom=0; jatom < this-> natom; ++jatom) {
            if (jatom == iatom) continue;
            xt::xtensor<double, 1> muij = (xt::col(dist, iatom) - xt::col(dist, jatom)) / xt::norm_l2(xt::row(xyz, iatom) - xt::row(xyz, jatom));
            xt::xtensor<double, 1> nuij{muij};
            xt::xtensor<double, 1>  sij{muij};

            if (this->biased)
            {
                double chi = BraggSlaterRadii.at(this->symbols[iatom]) / BraggSlaterRadii.at(this->symbols[jatom]);
                double uij = (chi - 1.) / (chi + 1.);
                double aij = uij / (uij * uij - 1.);
                if (aij > 0.5) aij = 0.5;
                else if (aij < -0.5) aij = -0.5;
                nuij = muij + aij * (1 - xt::square(muij));
            }

            beckegrid::sk(this->k, nuij, sij);
            s *= sij;
        }
        xt::row(pres, iatom) = s;
    }

    pres /= xt::sum(pres, {0});
    return pres;
}

xt::xtensor<double, 1> beckegrid::BeckeFuzzyCell::get_weight(const xt::xtensor<double, 2>& grid, const int iatom)
{
    auto pres = this->get_weight_p(grid);
    return xt::row(pres, iatom);
}

[[deprecated("this method returns one unnecessary array, use build_grid2 instead\n")]]
std::tuple<xt::xtensor<double, 4>, xt::xtensor<double, 3>> beckegrid::BeckeFuzzyCell::build_grid()
{
    this->xwleb = lebedev::lebedev_rule(this->nleb);                     // (4, nang)
    this->zcheb = xt::arange(1., static_cast<double>(this->nrad)+1.);   // (nrad, )
    this->xrwcheb = std::vector<xt::xtensor<double, 2>>{this->natom, xt::zeros<double>(xt::xtensor<double, 2>::shape_type{3, this->nrad})}; // [natom, (3, nrad)]
    size_t    nang = xwleb.shape(1);
    size_t natgrid = this->nrad * nang;
    xt::xtensor<double, 4>        grid = xt::zeros<double>({this->natom, this->nrad,  nang, 3+this->natom});            // (natom, nrad, nang, 3+natom) to be returned
    xt::xtensor<double, 3> grid_global = xt::zeros<double>(xt::xtensor<double, 3>::shape_type{this->natom, natgrid, 3}); // (natom, nradxnang, 3)        to be returned

    for (auto iatom=0; iatom < this->natom; ++iatom)
    {
        xt::xtensor<double, 2> zxrw = gaussCheby2(this->nrad, this->rms[iatom]);
        xt::row(this->xrwcheb[iatom], 0) = xt::row(zxrw, 1); // x
        xt::row(this->xrwcheb[iatom], 1) = xt::row(zxrw, 2); // r
        xt::row(this->xrwcheb[iatom], 2) = xt::row(zxrw, 3); // w

        for (auto j=0; j < this->nrad; ++j)
        {
            auto view1 = xt::view(grid, iatom, j, xt::all(), xt::range(0, 3));             // 坐标
            auto view2 = xt::view(grid, iatom, j, xt::all(), xt::range(3, 3+this->natom)); // 权重
            view1 = xt::transpose(xt::view(xwleb, xt::range(0,3), xt::all()) * zxrw(2, j) + xt::view(this->xyz, iatom, xt::all(), xt::newaxis()));
            view2 = xt::transpose(this->get_weight_p(view1));
        }

        this->weights[iatom] = xt::ravel<xt::layout_type::row_major>(xt::view(grid, iatom, xt::all(), xt::all(), 3+iatom));        // w
        xt::view(grid_global, iatom, xt::all(), 0) = xt::ravel<xt::layout_type::row_major>(xt::view(grid, iatom, xt::all(), xt::all(), 0)); // x
        xt::view(grid_global, iatom, xt::all(), 1) = xt::ravel<xt::layout_type::row_major>(xt::view(grid, iatom, xt::all(), xt::all(), 1)); // y
        xt::view(grid_global, iatom, xt::all(), 2) = xt::ravel<xt::layout_type::row_major>(xt::view(grid, iatom, xt::all(), xt::all(), 2)); // z
    }

    return {grid, grid_global};
}

std::vector<xt::xtensor<double, 2>> beckegrid::BeckeFuzzyCell::build_grid2()
{
    this->xwleb = lebedev::lebedev_rule(this->nleb); // (3+1, nang) lebedev
    const size_t    nang = xwleb.shape(1);
    const size_t natgrid = this->nrad * nang;
    this->zcheb = xt::arange(1., static_cast<double>(this->nrad)+1.); // (nrad, ) 所有原子共用一套
    this->xrwcheb = std::vector<xt::xtensor<double, 2>>{this->natom, xt::zeros<double>(xt::xtensor<double, 2>::shape_type{3, this->nrad})}; // [natom, (3, nrad)]
    this->weights = std::vector<xt::xtensor<double, 1>>{this->natom, xt::zeros<double>(xt::xtensor<double, 1>::shape_type{natgrid})}; // [natom, natgrid]
    
    std::vector<xt::xtensor<double, 2>> grid_global{this->natom, xt::zeros<double>(xt::xtensor<double, 2>::shape_type{natgrid, 3})}; // [natom, (nradxnang, 3)] to be returned
    
#pragma omp parallel for
    for (auto iatom=0; iatom < this->natom; ++iatom)
    {
        xt::xtensor<double, 2> zxrw = gaussCheby2(this->nrad, this->rms[iatom]); // (4, nrad)
        xt::row(this->xrwcheb[iatom], 0) = xt::row(zxrw, 1); // x
        xt::row(this->xrwcheb[iatom], 1) = xt::row(zxrw, 2); // r
        xt::row(this->xrwcheb[iatom], 2) = xt::row(zxrw, 3); // w

        for (auto j=0; j < this->nrad; ++j)
        {
            auto view_tmp = xt::view(grid_global[iatom], xt::range(j * nang, (j+1) * nang), xt::all()); // view shape (natgrid, 3)
            view_tmp = xt::transpose(zxrw(2, j) * xt::view(xwleb, xt::range(0, 3), xt::all()) + xt::view(this->xyz, iatom, xt::all(), xt::newaxis())); // assign atgrid coord
        }

        this->weights[iatom] = this->get_weight(grid_global[iatom], iatom);
    }

    return grid_global;
}

