#include "cspline.hpp"


cspline::CubicSpline::CubicSpline(const xt::xtensor<double, 1>& xx, const xt::xtensor<double, 1>& yy)
{
    assert(xx.size() == yy.size());
    this->x = xx;
    this->y = yy;
    this->spline_generated = false;
}

void cspline::CubicSpline::generate_spline()
{
    this->spline = std::vector<std::vector<double>>(this->x.size()-1,std::vector<double>(4,0.0));
    std::vector<double> A(this->x.size(),0.0);
    std::vector<double> B(this->x.size(),0.0);
    std::vector<double> C(this->x.size(),0.0);
    std::vector<double> D(this->x.size(),0.0);
    std::vector<double> Y(this->x.size(),0.0);

    double h0, h1, r0, r1;
    h0   =  this->x[1] - this->x[0];
    h1   =  this->x[2] - this->x[1];
    r0   = (this->y[1] - this->y[0]) / h0;
    r1   = (this->y[2] - this->y[1]) / h1;
    B[0] = h1 * (h0 + h1);
    C[0] = (h0 + h1) * (h0 + h1);
    Y[0] = r0 * (3 * h0 * h1 + 2 * h1 * h1) + r1 * h0 * h0;

    for(auto i=1; i < this->x.size()-1; ++i) {
        h0   =  this->x[i]   - this->x[i-1];
        h1   =  this->x[i+1] - this->x[i];
        r0   = (this->y[i]   - this->y[i-1]) / h0;
        r1   = (this->y[i+1] - this->y[i])   / h1;
        A[i] = h1;
        B[i] = 2 * (h0 + h1);
        C[i] = h0;
        Y[i] =3 * (r0 * h1 + r1 * h0);
    }

    A[this->x.size()-1] = (h0 + h1) * (h0 + h1);
    B[this->x.size()-1] = h0 * (h0 + h1);
    Y[this->x.size()-1] = r0 * h1 * h1 + r1 * (3 * h0 * h1 + 2 * h0 * h0);

    C[0] = C[0] / B[0];
    for(auto i=1; i < this->x.size()-1; ++i) C[i] = C[i] / (B[i] - A[i] * C[i-1]);

    Y[0] = Y[0] / B[0];
    for(auto i=1; i < this->x.size(); ++i) Y[i] = (Y[i] - A[i] * Y[i-1]) / (B[i] - A[i] * C[i-1]);

    D[this->x.size()-1] = Y[this->x.size()-1];
    for(auto i=this->x.size()-1; i > 0; i--) D[i-1] = Y[i-1] - C[i-1] * D[i];

    double dx, dy;
    for(auto i=0; i < this->x.size()-1; ++i) {
        dx                 = 1. / (this->x[i+1] - this->x[i]);
        dy                 = (this->y[i+1] - this->y[i]) * dx;
        this->spline[i][0] = this->y[i]; //ai
        this->spline[i][1] = D[i]; //bi
        this->spline[i][2] = dx * (3 * dy - 2 * D[i] - D[i+1]); //ci
        this->spline[i][3] = dx * dx * (-2 * dy + D[i] + D[i+1]); //di
    }

    this->spline_generated = true;
}

double cspline::CubicSpline::eval(const double xx) const
{
    assert(this->spline_generated);
    if (xx < this->x[0] || xx > this->x[this->x.size()-1]) {std::cerr << "Out of bounds!" << std::endl; exit(-1);}
    // 没有外推
    // if(xx < this->x[0]) return this->y.front();
    // if(xx >= this->x.back()) return this->y.back();

    for(auto i=1; i < this->x.size(); ++i)
    {
        if(xx <= this->x[i])
        {
            double xpar = (xx - this->x[i-1]);
            return this->spline[i-1][0] + this->spline[i-1][1] * xpar + this->spline[i-1][2] * xpar * xpar + this->spline[i-1][3] * xpar * xpar * xpar;
        }
    }

    return 0.;
}

