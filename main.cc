#include <cstdlib> // to get env
#include <iostream>
#include <cstring>
#include <cassert>

#ifdef TIMEIT
#include <chrono>
#endif

#include "include/dft.hpp"
#include "include/hf.hpp"


void test_hf()
{
    const char* datadir = std::getenv("DATADIR");
    if (datadir == nullptr) {std::cerr << "DATADIR not defined\n"; exit(-1);}
    const std::string LDATADIR{datadir};

    std::cout << LDATADIR << std::endl;

    libint2::Shell::do_enforce_unit_normalization(false);
    libint2::initialize();

    std::string xyzfile = std::string{LDATADIR} + "/xyz/ch3oh.xyz";
    std::string name    = std::string{LDATADIR} + "/gbs/cc-pvdz.g94";
    sf::HF::HF hf{xyzfile, name};

    std::cout << "\n!!! Running HF in Spherical form !!!\n";
    hf.scf(30, 1e-8, 1e-8, -1, "core", true);

    libint2::finalize();
}


void test_dft(std::string& mol)
{
    const char* datadir = std::getenv("DATADIR");
    if (datadir == nullptr) {std::cerr << "DATADIR not defined\n"; exit(-1);}
    const std::string LDATADIR{datadir};

    std::cout << LDATADIR << std::endl;

    libint2::Shell::do_enforce_unit_normalization(false);
    libint2::initialize();

    // std::string xyzfile = std::string{LDATADIR} + "/xyz/ch4.xyz";
    std::string xyzfile = mol;
    std::string name    = std::string{LDATADIR} + "/gbs/cc-pvdz.g94";
    sf::DFT::DFT dft{xyzfile, name};

    dft.init();


    std::cout << "\n!!! Running DFT in Spherical form !!!\n";
    dft.scf(30, 1e-8, 1e-8, -1, "core", 1, 8, true);

    libint2::finalize();
}



int main(int argc, char** argv)
{
    assert(argc == 2);
    std::string mol{argv[1]};

    #ifdef TIMEIT
    const auto t0{std::chrono::steady_clock::now()};
    #endif
    test_dft(mol);
    // test_hf();
    #ifdef TIMEIT
    const auto t1{std::chrono::steady_clock::now()};
    std::cout << "\nNormal termination or failed termination, you spent " << std::chrono::duration<double>{t1-t0} << " sec on it." << std::endl;
    #endif
}
