**[中文](README_zh.md)**

# SegFault

## About

Low-performance, hard-to-use quantum chemistry C++ code, currently in its very early stage and only supports closed-shell HF, LDA and GGA calculation. 

The code is designed to be input-file-free (or use .json file to pass tasks in future), users should use this code as third-party library and incorporate it into computational tasks and use corresponding APIs. 

## Dependencies

1. [xtensor](https://xtensor.readthedocs.io/en/latest/) (header only)
2. [libint-2.11.2](https://github.com/evaleev/libint)
3. [libxc-7.0.0](https://gitlab.com/libxc/libxc)
4. (optional but strongly recommended) intel mkl and intel icpx compiler (version 2025.3).

## Installation
```bash
source /path/to/intel/oneapi/setvars.sh
./worker.sh
cd install/bin/
./segfault /path/to/mol.xyz /path/to/gbs [X_id C_id]
```

## Example
see main.cc

## Contributing

This project welcomes (and urgently needs) contributions of any kinds, include but not not limited to adding new features, performance optimization, bug fixing, docs update, etc. 
