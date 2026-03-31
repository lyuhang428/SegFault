#include <iostream>
#include <format>
#include <cassert>
#include <cstring>
#include <map>
#include <array>

#ifdef TIMEIT
#include <chrono>
#endif

#include "include/dft.hpp"


static const std::map<std::string, std::array<int, 2>> xc = {
    {"pz", {1, 9}},
    {"svwn5", {1, 8}}, 
    {"pbe", {101, 130}}, 
    {"blyp", {106, 131}},
    {"pbesol", {116, 133}},
    {"pw91", {109, 134}}
};




void test_dft(std::string xyzfile, std::string gbsfile, int X_id, int C_id)
{
    libint2::Shell::do_enforce_unit_normalization(false);
    libint2::initialize();

    sf::DFT::DFT dft{xyzfile, gbsfile};
    dft.scf(30, 1e-8, 1e-8, -1, "core", X_id, C_id, true, 75, 29, 3, true, "becke");

    libint2::finalize();
}



int main(int argc, char** argv)
{
    assert(argc == 5);
    std::string xyzfile{argv[1]};
    std::string gbsfile{argv[2]};
    int X_id = std::stoi(argv[3]);
    int C_id = std::stoi(argv[4]);

    std::cout << xyzfile << std::endl;
    std::cout << gbsfile << std::endl;

#ifdef TIMEIT
const auto t0{std::chrono::steady_clock::now()};
#endif

    test_dft(xyzfile, gbsfile, X_id, C_id);

#ifdef TIMEIT
const auto t1{std::chrono::steady_clock::now()};
const auto wtime{std::chrono::duration<float>(t1-t0)};
std::cout << "wtime " << wtime << std::endl;
#endif


    std::cout << "\nsafe here\n" << std::endl;
}
