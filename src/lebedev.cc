#include "../include/lebedev.hpp"

xt::xtensor<double, 2> lebedev::lebedev_rule(const int n)
{
    auto it = std::find(std::begin(LEBEDEV_ORDER), std::end(LEBEDEV_ORDER), n);
    if (it == std::end(LEBEDEV_ORDER)) {std::cerr << "Error: Requested Lebedev order not possible" << std::endl; exit(-1);}
    
    int index = it - std::begin(LEBEDEV_ORDER);
    int  nleb = LEBEDEV_LEVEL[index];
    xt::xtensor<double, 2> xwleb{{4, static_cast<size_t>(nleb)}, 0.}; // to be returned
    double theta, phi, w;

    std::stringstream iss;
    iss << std::setw(3) << std::setfill('0') << n;
    std::string filename = "../../data/lebedev/lebedev_" + iss.str() + ".dat";
    std::ifstream ifile(filename, std::ios::in);
    if (!ifile.is_open()) {std::cerr << "Error: File open failed" << std::endl; exit(-1);}
    std::string buffer;
    
    for (size_t id=0; id < nleb; ++id)
    {
        if (!std::getline(ifile, buffer)) {std::cerr << "Error: failed to read file" << std::endl; exit(-1);}
        std::istringstream isss(buffer);
        isss >> phi >> theta >> w;
        theta = theta / 180. * M_PI;
        phi   = phi   / 180. * M_PI;
        xwleb(0, id) = std::sin(theta) * std::cos(phi);
        xwleb(1, id) = std::sin(theta) * std::sin(phi);
        xwleb(2, id) = std::cos(theta);
        xwleb(3, id) = w * 4. * M_PI;
    }
    return xwleb;
}


