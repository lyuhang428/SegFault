**[中文](README_zh.md)**

# SegFault

## About

Low-performance, hard-to-use quantum chemistry C++ code, currently in its very early stage and only supports closed-shell Hartree-Fock and LDA level DFT calculation. 

The code is designed to be input-file-free (or use .json file to pass tasks), users should use this code as third-party library and incorporate it into computational tasks and use corresponding APIs. 

## Dependencies

1. [Eigen](https://libeigen.gitlab.io/eigen/docs-5.0/) (header only)
2. [xtensor](https://xtensor.readthedocs.io/en/latest/) (header only)
3. [libint2](https://github.com/evaleev/libint)
4. [libxc](https://gitlab.com/libxc/libxc)
5. [gsl](https://www.gnu.org/software/gsl/)
6. (optional) intel mkl and intel icpx compiler

## Installation
```bash
source /path/to/intel/oneapi/setvars.sh
export DATADIR=path/to/data/
cmake -B build -DCMAKE_CXX_COMPILER=icpx -DCMAKE_INSTALL_PREFIX=./install
cmake --build build/ -j8
cd build && make install
cd ../install/bin/
./segfault
```

## Example
see main.cc

## Contributing

This project welcomes (and urgently needs) contributions of any kind, include but not not limited to adding new features, performance optimization, bug fixing, docs update, etc. 
