#include <cstdlib>
#include <iostream>

#include "include/hf.hpp"
#include "include/dft.hpp"


void test_hf()
{
    const char* datadir = std::getenv("DATADIR");
    if (datadir == nullptr) {std::cerr << "DATADIR not defined\n"; exit(-1);}
    const std::string LDATADIR{datadir};
    
    std::cout << LDATADIR << std::endl;

    libint2::Shell::do_enforce_unit_normalization(false);
    libint2::initialize();

    std::string xyzfile = std::string{LDATADIR} + "/xyz/co2.xyz";
    std::string name    = std::string{LDATADIR} + "/gbs/cc-pvdz.g94";
    sf::HF::HF hf{xyzfile, name};

    std::cout << "\n!!! Running HF in Spherical form !!!\n";
    hf.scf(30, 1e-6, 1e-6, -1, "core", true);
    

    libint2::finalize();
}

void test_dft()
{
    const char* datadir = std::getenv("DATADIR");
    if (datadir == nullptr) {std::cerr << "DATADIR not defined\n"; exit(-1);}
    const std::string LDATADIR{datadir};

    std::cout << LDATADIR << std::endl;

    libint2::Shell::do_enforce_unit_normalization(false);
    libint2::initialize();

    std::string xyzfile = std::string{LDATADIR} + "/xyz/ch3oh.xyz";
    std::string name    = std::string{LDATADIR} + "/gbs/cc-pvdz.g94";
    sf::DFT::DFT dft{xyzfile, name};

    dft.init();


    std::cout << "\n!!! Running DFT in Spherical form !!!\n";
    dft.scf(30, 1e-8, 1e-8, -1, "core", 1, 8, true);

    libint2::finalize();
}



int main()
{
    
    test_hf();
    std::cout << std::string(60, '-') << std::endl;
    test_dft();
}
