<h1 id="english">SegFault</h1>
<a href="#chinese">[中文]</a>

## About

Independent electronic structure code development (GTO basis). This project is currently in its very early stage and only supports closed-shell Hartree-Fock, LDA, GGA, and meta-GGA calculations. 

The code is designed to be *input-file-free* (or interface to Fortran `namelist` to pass tasks in future). Users should use this code as a library and incorporate it into their own computational tasks and call corresponding APIs. 

## Dependencies

1. [xtensor](https://xtensor.readthedocs.io/en/latest/) (header only)
2. [libint-2.11.2](https://github.com/evaleev/libint)
3. [libxc-7.0.0](https://gitlab.com/libxc/libxc)
4. `Intel MKL` and intel icpx compiler (version 2025.3). Interface to `GNU` compiler and `OpenBLAS` in progress.


## Code structure
```txt
include/
├── beckegrid.hpp
├── constants.hpp
├── dft.hpp
├── hf.hpp
├── lblas.hpp
├── lebedev.hpp
└── mints.hpp

src/
├── beckegrid.cc
├── dft.cc
├── hf.cc
├── lebedev.cc
└── mints.cc

data/
├── gbs
├── lebedev
└── xyz
```

- `lblas.hpp` interfaces to `Intel MKL` to perform linear algebra.
- `beckgrid.cc` build Becke grid and grid weight for exchange-correlation numerical quadrature.
- `lebedev.cc` use Lebedev quadrature scheme to do numerical quadrature over solid angle.
- `mints.cc` interface to `libint` to deal with molecular integral.
- `dft.cc` contain the main routines to perform DFT calculations.
- `hf.cc` contain the main routines to perform Hartree-Fock calculations.
- `data/gbs` contain basis set in Gaussian format.
- `data/xyz` contain molecules xyz file.
- `data/levedev` contain the roots and weights of Levedev quadrature.

## Installation
```bash
source /path/to/intel/oneapi/setvars.sh
./worker.sh
cd ./install/bin/
./segfault /path/to/mol.xyz /path/to/gbs [X_id C_id] # X_id is exchange functional id in `libxc`; C_id is correlation functional id in `libxc`
```

## Example
```c++
// perform DFT calculations
#include "include/dft.hpp"
libint2::Shell::do_enforce_unit_normalization(false); // initialize libint
libint2::initialize();
sf::DFT::DFT dft{xyzfile, gbsfile};
dft.scf(50,       // maxiter
        1e-8,     // e_convergence
        1e-8,     // d_convergence
        -1,       // nbuffer for DIIS
        "core",   // initial guess
        X_id,     // exchange functional id in libxc
        C_id,     // correlation functional id in libxc
        true,     // use spherical basis function
        75,       // number of radial grid
        29,       // Lebedev quadrature order
        3,        // 3-order switching function
        true,     // biased Becke partition
        "becke"); // finite to quasi-infinite mapping scheme

libint2::finalize(); // close libint

// perform HF calculations
#include "include/hf.hpp"
libint2::Shell::do_enforce_unit_normalization(false); // initialize libint
libint2::initialize();
sf::HF::HF hf{xyzfile, gbsfile};
hf.scf(50,     // maxiter
       1e-8,   // e_convergence
       1e-8,   // d_convergence
       -1,     // nbuffer for DIIS
       "core", // initial guess
       true);  // use spherical basis function

libint2::finalize(); // close libint
```

## Validation/Benchmark
Validation against Psi4 (version 1.9.1)
Machine: ThinkPad T14p, 13th Gen Intel(R) Core(TM) i5-13500H, 16 threads.
cc-pvdz basis set.
- **PBE**

| System | $\quad\Delta E_\text{tot}\quad$ | $\quad\Delta E_\text{1e}\quad$ | $\quad\Delta E_\text{2e}\quad$ | $\quad\Delta E_\text{xc}\quad$ | $\Delta E_\text{HOMO-LUMO gap}$ | Wall time (sec) | 
| :- | :-: | :-: | :-: | :-: | :-: | -: |
| methylamine | -3.2703e-06 | -2.5987e-06 | 2.3825e-06 | -3.0727e-06 | 2.6090e-07 | 6.1806 |
| carbon_monoxide | 6.5880e-07 | 1.3818e-05 | -1.5620e-05 | 2.4534e-06 | 1.9010e-07 | 1.4942 |
| benzene | 6.7156e-06 | 8.3642e-05 | -8.9751e-05 | 1.2736e-05 | -6.6730e-07 | 41.8788 |
| 1-hexanol | -4.1803e-05 | -1.2538e-04 | 1.3985e-04 | -5.6411e-05 | -1.5520e-06 | 162.4541 |
| oxirane | -5.0822e-07 | -3.9450e-05 | 4.3488e-05 | -4.5791e-06 | 9.7150e-07 | 9.8304 |
| pyrrole | -4.1548e-06 | 1.0334e-05 | -1.2756e-05 | -1.8040e-06 | 5.5700e-08 | 26.9344 |
| furan | -3.0730e-06 | 1.5352e-05 | -1.5990e-05 | -2.5054e-06 | 4.3170e-07 | 27.8473 |
| propane | -1.2031e-05 | -6.1754e-05 | 6.7415e-05 | -1.7728e-05 | 5.7680e-07 | 18.5657 |
| 1-propanol | 1.0036e-05 | 1.9686e-05 | -2.4752e-05 | 1.5043e-05 | -1.5490e-06 | 67.8929 |
| 1-butylamine | 3.0350e-07 | -2.6457e-04 | 2.8795e-04 | -2.3160e-05 | -3.0931e-06 | 66.9922 |
| water | 1.1774e-07 | 5.2071e-05 | -5.9104e-05 | 7.1467e-06 | -3.3560e-07 | 1.0797 |
| 1-hexylamine | -6.8038e-06 | -1.0323e-04 | 1.1422e-04 | -1.7932e-05 | 1.1438e-06 | 176.0360 |
| methanol | -1.7309e-06 | -3.1290e-06 | 3.4408e-06 | -2.0605e-06 | 6.4805e-06 | 4.1867 |
| ethene | -5.7397e-07 | 4.1682e-05 | -4.6175e-05 | 3.9049e-06 | 1.8240e-07 | 3.3631 |
| methane | 2.5309e-08 | -4.5112e-06 | 5.3238e-06 | -7.9313e-07 | 9.0110e-07 | 1.9833 |
| ethylamine | -7.6544e-06 | -4.9935e-05 | 5.5777e-05 | -1.3533e-05 | -9.1540e-07 | 18.0129 |
| ethanol | 1.1160e-05 | -4.8624e-05 | 5.6853e-05 | 2.8938e-06 | 7.5760e-06 | 14.6105 |
| hexane | 1.4781e-05 | 2.8967e-05 | -3.1717e-05 | 1.7421e-05 | -1.1650e-07 | 130.9111 |
| cyclopentane | 1.9080e-05 | 8.7994e-05 | -9.3962e-05 | 2.4965e-05 | 2.6800e-08 | 54.5885 |
| butane | 1.2717e-05 | -3.5731e-05 | 3.7339e-05 | 1.1050e-05 | 5.9610e-07 | 41.5918 |
| 1-pentanol | -5.7729e-05 | -9.9320e-05 | 1.2619e-04 | -8.4705e-05 | 5.0051e-06 | 102.8831 |
| 1-butanol | 2.6880e-05 | -9.4562e-05 | 1.0508e-04 | 1.6285e-05 | 9.2200e-07 | 74.2747 |
| cyclopropane | 1.1728e-05 | 4.1420e-05 | -4.5669e-05 | 1.5943e-05 | -9.2150e-07 | 11.6828 |
| hydrogen_peroxide | -2.2752e-07 | 2.0495e-05 | -2.3157e-05 | 2.4183e-06 | 2.3700e-07 | 2.3459 |
| ammonia | 3.7417e-07 | 2.9976e-05 | -3.3210e-05 | 3.6035e-06 | 9.6700e-08 | 1.6706 |
| cyclobutane | -1.0257e-05 | 1.7592e-04 | -2.0002e-04 | 1.3786e-05 | -1.0660e-05 | 22.9783 |
| pentane | -5.5838e-06 | -2.1463e-04 | 2.2873e-04 | -1.9761e-05 | 2.0828e-06 | 79.4211 |
| pyridine | 4.5146e-06 | -4.4350e-05 | 4.7036e-05 | 1.7380e-06 | 2.0220e-07 | 42.1721 |
| 1-propylamine | 1.1753e-05 | 6.6638e-05 | -7.4420e-05 | 1.9476e-05 | -3.9660e-07 | 51.0550 |
| 1-propene | 2.9772e-07 | 7.4772e-06 | -8.6289e-06 | 1.4183e-06 | -7.0850e-07 | 13.0076 |
| carbon_dioxide | 1.4071e-06 | -3.5163e-05 | 4.2448e-05 | -5.8988e-06 | -8.6200e-08 | 1.9148 |
| cyclohexane | -4.9942e-05 | 2.8359e-04 | -3.3186e-04 | -1.7792e-06 | -3.0933e-05 | 84.4095 |
| 1-pentene | -3.4967e-05 | -2.7088e-05 | 3.3169e-05 | -4.1124e-05 | -3.2000e-08 | 64.0157 |
| 1-pentylamine | -9.2611e-06 | -6.9553e-05 | 7.6401e-05 | -1.6219e-05 | -1.9798e-06 | 129.7257 |
| 1-butene | 1.1407e-05 | -1.2461e-04 | 1.5473e-04 | -1.8756e-05 | 1.4308e-06 | 38.8057 |
| ethane | 2.5769e-06 | -1.5856e-05 | 1.7729e-05 | 6.8581e-07 | -6.0010e-07 | 4.9346 |
| 1-hexene | 2.8032e-06 | -1.2182e-05 | 8.5036e-06 | 6.3801e-06 | 1.4920e-07 | 136.3154 |

- **BLYP**

| System | $\quad\Delta E_\text{tot}\quad$ | $\quad\Delta E_\text{1e}\quad$ | $\quad\Delta E_\text{2e}\quad$ | $\quad\Delta E_\text{xc}\quad$ | $\Delta E_\text{HOMO-LUMO gap}$ | Wall time (sec) | 
| :- | :-: | :-: | :-: | :-: | :-: | -: |
| methylamine | -3.4158e-06 | -5.4660e-06 | 5.3908e-06 | -3.3592e-06 | 2.2086e-06 | 8.1370 |
| carbon_monoxide | 7.6566e-07 | -1.1822e-04 | 1.3371e-04 | -1.4737e-05 | -2.4410e-07 | 2.8111 |
| benzene | 5.5895e-06 | 7.0633e-05 | -7.6796e-05 | 1.1663e-05 | 2.0710e-07 | 41.9015 |
| 1-hexanol | -3.7437e-05 | -1.7841e-05 | 2.7389e-05 | -4.7127e-05 | -7.6600e-07 | 156.7170 |
| oxirane | -9.4525e-07 | -7.0850e-06 | 7.8288e-06 | -1.7219e-06 | 5.9990e-07 | 12.8891 |
| pyrrole | -3.6845e-06 | 2.5182e-05 | -2.8764e-05 | -1.7326e-07 | 2.3940e-07 | 31.0708 |
| furan | -2.4140e-06 | 1.3709e-05 | -1.4085e-05 | -2.1086e-06 | -1.4460e-07 | 30.7187 |
| propane | -1.1635e-05 | -1.1520e-04 | 1.2569e-04 | -2.2161e-05 | -6.2900e-07 | 21.5090 |
| 1-propanol | 1.1320e-05 | 4.9515e-06 | -7.5846e-06 | 1.3894e-05 | -1.4827e-06 | 79.1883 |
| 1-butylamine | 4.9500e-07 | 2.9983e-05 | -3.2166e-05 | 2.5956e-06 | -2.3196e-06 | 77.7563 |
| water | 2.1403e-07 | -1.5819e-04 | 1.7886e-04 | -2.0468e-05 | 3.4300e-08 | 1.7797 |
| 1-hexylamine | -4.4126e-06 | -1.5205e-04 | 1.6478e-04 | -1.7279e-05 | 5.2580e-07 | 176.4170 |
| methanol | -1.7133e-06 | -6.6279e-06 | 8.2219e-06 | -3.3250e-06 | 4.2416e-06 | 6.9046 |
| ethene | -2.8329e-07 | 1.4407e-04 | -1.6367e-04 | 1.9302e-05 | -5.8570e-07 | 5.2840 |
| methane | 4.3539e-08 | -9.4489e-06 | 1.0604e-05 | -1.1180e-06 | 2.7670e-07 | 3.1808 |
| ethylamine | -7.4635e-06 | -7.0719e-05 | 7.9068e-05 | -1.5850e-05 | -4.7390e-07 | 21.2482 |
| ethanol | 1.3092e-05 | -6.8260e-05 | 7.7500e-05 | 3.8160e-06 | 6.7762e-06 | 18.5641 |
| hexane | 1.3382e-05 | 5.1310e-05 | -5.4904e-05 | 1.6867e-05 | -1.4446e-06 | 134.1697 |
| cyclopentane | 1.6050e-05 | 8.8102e-05 | -9.4682e-05 | 2.2547e-05 | 8.3010e-07 | 57.1636 |
| butane | 9.7959e-06 | 2.2874e-05 | -2.5640e-05 | 1.2503e-05 | 3.6290e-07 | 44.3323 |
| 1-pentanol | -4.8670e-05 | 1.3929e-05 | 2.9612e-06 | -6.5672e-05 | 3.2538e-06 | 111.4073 |
| 1-butanol | 2.7061e-05 | 6.8702e-05 | -6.9582e-05 | 2.7858e-05 | -9.7700e-07 | 84.2336 |
| cyclopropane | 9.2693e-06 | 4.0419e-05 | -4.4120e-05 | 1.2937e-05 | -3.2500e-07 | 14.8105 |
| hydrogen_peroxide | -4.3466e-07 | 7.9258e-05 | -8.9659e-05 | 9.9500e-06 | 8.2650e-07 | 4.2311 |
| ammonia | 4.0432e-07 | 1.5176e-05 | -1.6666e-05 | 1.8891e-06 | 3.6860e-07 | 2.7723 |
| cyclobutane | -1.1024e-05 | 1.3205e-04 | -1.5119e-04 | 8.0648e-06 | -3.4779e-06 | 28.4853 |
| pentane | -4.1732e-06 | -1.4490e-04 | 1.5469e-04 | -1.4052e-05 | 1.3412e-06 | 82.2843 |
| pyridine | 2.3877e-06 | -6.2271e-05 | 6.6459e-05 | -1.8908e-06 | 7.8800e-07 | 44.7745 |
| 1-propylamine | 1.0828e-05 | -9.5405e-05 | 9.8680e-05 | 7.4945e-06 | 1.3623e-06 | 57.1780 |
| 1-propene | -6.0543e-07 | 5.1962e-06 | -6.4260e-06 | 5.9320e-07 | -8.9520e-07 | 17.6234 |
| carbon_dioxide | 1.6545e-06 | -5.6646e-05 | 6.8060e-05 | -9.7816e-06 | 1.1220e-07 | 4.6161 |
| cyclohexane | -2.7313e-05 | 2.2525e-04 | -2.6977e-04 | 1.7092e-05 | -2.7817e-05 | 87.1533 |
| 1-pentene | -2.9387e-05 | 1.4313e-05 | -1.1801e-05 | -3.1975e-05 | -2.6790e-07 | 65.2221 |
| 1-pentylamine | -8.7238e-06 | 1.9610e-05 | -1.8562e-05 | -9.8809e-06 | -2.6528e-06 | 127.4520 |
| 1-butene | 1.4674e-05 | -1.1320e-04 | 1.4084e-04 | -1.3023e-05 | 2.1204e-06 | 42.7310 |
| ethane | 2.1856e-06 | -2.2522e-05 | 2.5187e-05 | -4.9816e-07 | -4.8170e-07 | 8.8190 |
| 1-hexene | 3.2542e-06 | 3.2270e-05 | -3.3780e-05 | 4.6620e-06 | -8.2770e-07 | 142.2336 |

- **HF**

| System | $\quad\Delta E_\text{tot}\quad$ | Wall time (sec) | 
| :- | :-: | -: |
| methanol | 3.1604e-10 | 3.003016 | 
| ethanol | 3.5021e-10 | 13.171072 | 
| 1-propanol | 6.8701e-10 | 35.912599 | 
| 1-butanol | 2.6955e-10 | 96.356437 | 
| 1-pentanol | -5.6957e-11 | 195.1718 | 
| 1-hexanol | 2.0833e-10 | 353.214047 | 
| methane | 4.5851e-11 | 0.945244 | 
| ethane | -9.1902e-11 | 5.400898 | 
| propane | -5.9757e-11 | 18.332309 | 
| butane | 2.2504e-10 | 5.672531 | 
| pentane | -5.0156e-10 | 124.540576 | 
| hexane | -4.2908e-10 | 252.264091 | 
| ethene | 2.9544e-10 | 3.192671 | 
| 1-propene | -1.7594e-10 | 11.814476 | 
| 1-butene | -5.4098e-10 | 39.080278 | 
| 1-pentene | -6.9781e-10 | 92.765139 | 
| 1-hexene | -9.8126e-10 | 187.272921 | 
| methylamine | 1.3522e-10 | 4.016299 | 
| ethylamine | 9.5213e-11 | 16.972886 | 
| 1-propylamine | 2.1083e-10 | 48.243879 | 
| 1-butylamine | 7.1347e-10 | 110.805998 | 
| 1-pentylamine | 2.4392e-10 | 219.96451 | 
| 1-hexylamine | -4.2195e-10 | 389.22106 | 
| cyclopropane | 8.4412e-11 | 12.852455 | 
| cyclobutane | 5.5962e-11 | 38.54993 | 
| cyclopentane | -1.4057e-10 | 91.709966 | 
| cyclohexane | 1.3378e-10 | 192.254076 | 
| water | 3.1559e-10 | 0.438969 | 
| ammonia | 3.7608e-10 | 0.632635 | 
| pyrrole | 5.8455e-10 | 40.915752 | 
| pyridine | 8.0149e-11 | 72.489318 | 
| oxirane | 8.2983e-10 | 7.919546 | 
| hydrogen_peroxide | -6.6203e-10 | 1.643334 | 
| furan | 1.1821e-09 | 35.161369 | 
| carbon_monoxide | 3.5479e-09 | 0.761816 | 
| carbon_dioxide | 5.7953e-09 | 2.485237 | 


## Contributing

This project welcomes contributions of any kinds, include but not not limited to adding new features, performance optimization, bug fixing, docs update, etc. 

## Future development
Implement routines for force calculation.

---

<h1 id="chinese">SegFault</h1>
<a href="#english">[English]</a>

## 简介

电子结构代码开发项目（高斯基组）。该项目目前处于非常早期阶段，仅支持闭壳层 Hartree-Fock, LDA, GGA, 以及 meta-GGA 计算。 

设计逻辑是不使用输入文件（或未来接口 Fortran `namelist` 传递任务）。建议用户将代码作为第三方库集成到计算任务中并依据需求使用相应的接口。

## 依赖库

1. [xtensor](https://xtensor.readthedocs.io/en/latest/) (header only)
2. [libint-2.11.2](https://github.com/evaleev/libint)
3. [libxc-7.0.0](https://gitlab.com/libxc/libxc)
4. `Intel MKL` and intel icpx compiler (version 2025.3). Interface to `GNU` compiler and `OpenBLAS` in progress.

## 代码结构
```txt
include/
├── beckegrid.hpp
├── constants.hpp
├── dft.hpp
├── hf.hpp
├── lblas.hpp
├── lebedev.hpp
└── mints.hpp

src/
├── beckegrid.cc
├── dft.cc
├── hf.cc
├── lebedev.cc
└── mints.cc

data/
├── gbs
├── lebedev
└── xyz
```

- `lblas.hpp` 提供 `Intel MKL` 接口用于线性代数计算。
- `beckgrid.cc` 构造 Becke 网格以及网格权重用于交换相关势矩阵元计算。
- `lebedev.cc` Lebedev 球面求积。
- `mints.cc` 接口 `libint` 计算分子积分。
- `dft.cc` 进行 DFT 计算的主程序。
- `hf.cc` 进行 HF 计算的主程序。
- `data/gbs` 高斯基组。
- `data/xyz` 分子坐标。
- `data/levedev` Levedev 求积单位球网格坐标与权重。

## 安装
```bash
source /path/to/intel/oneapi/setvars.sh
./worker.sh
cd ./install/bin/
./segfault /path/to/mol.xyz /path/to/gbs [X_id C_id] # X_id is exchange functional id in `libxc`; C_id is correlation functional id in `libxc`
```

## 示例
```c++
// perform DFT calculations
#include "include/dft.hpp"
libint2::Shell::do_enforce_unit_normalization(false); // initialize libint
libint2::initialize();
sf::DFT::DFT dft{xyzfile, gbsfile};
dft.scf(50,       // maxiter
        1e-8,     // e_convergence
        1e-8,     // d_convergence
        -1,       // nbuffer for DIIS
        "core",   // initial guess
        X_id,     // exchange functional id in libxc
        C_id,     // correlation functional id in libxc
        true,     // use spherical basis function
        75,       // number of radial grid
        29,       // Lebedev quadrature order
        3,        // 3-order switching function
        true,     // biased Becke partition
        "becke"); // finite to quasi-infinite mapping scheme

libint2::finalize(); // close libint

// perform HF calculations
#include "include/hf.hpp"
libint2::Shell::do_enforce_unit_normalization(false); // initialize libint
libint2::initialize();
sf::HF::HF hf{xyzfile, gbsfile};
hf.scf(50,     // maxiter
       1e-8,   // e_convergence
       1e-8,   // d_convergence
       -1,     // nbuffer for DIIS
       "core", // initial guess
       true);  // use spherical basis function

libint2::finalize(); // close libint
```

## 结果验证
计算结果与 Psi4 (version 1.9.1) 进行比较。
硬件：ThinkPad T14p, 13th Gen Intel(R) Core(TM) i5-13500H, 16 线程.
基组：cc-pvdz
- **PBE**

| System | $\quad\Delta E_\text{tot}\quad$ | $\quad\Delta E_\text{1e}\quad$ | $\quad\Delta E_\text{2e}\quad$ | $\quad\Delta E_\text{xc}\quad$ | $\Delta E_\text{HOMO-LUMO gap}$ | Wall time (sec) | 
| :- | :-: | :-: | :-: | :-: | :-: | -: |
| methylamine | -3.2703e-06 | -2.5987e-06 | 2.3825e-06 | -3.0727e-06 | 2.6090e-07 | 6.1806 |
| carbon_monoxide | 6.5880e-07 | 1.3818e-05 | -1.5620e-05 | 2.4534e-06 | 1.9010e-07 | 1.4942 |
| benzene | 6.7156e-06 | 8.3642e-05 | -8.9751e-05 | 1.2736e-05 | -6.6730e-07 | 41.8788 |
| 1-hexanol | -4.1803e-05 | -1.2538e-04 | 1.3985e-04 | -5.6411e-05 | -1.5520e-06 | 162.4541 |
| oxirane | -5.0822e-07 | -3.9450e-05 | 4.3488e-05 | -4.5791e-06 | 9.7150e-07 | 9.8304 |
| pyrrole | -4.1548e-06 | 1.0334e-05 | -1.2756e-05 | -1.8040e-06 | 5.5700e-08 | 26.9344 |
| furan | -3.0730e-06 | 1.5352e-05 | -1.5990e-05 | -2.5054e-06 | 4.3170e-07 | 27.8473 |
| propane | -1.2031e-05 | -6.1754e-05 | 6.7415e-05 | -1.7728e-05 | 5.7680e-07 | 18.5657 |
| 1-propanol | 1.0036e-05 | 1.9686e-05 | -2.4752e-05 | 1.5043e-05 | -1.5490e-06 | 67.8929 |
| 1-butylamine | 3.0350e-07 | -2.6457e-04 | 2.8795e-04 | -2.3160e-05 | -3.0931e-06 | 66.9922 |
| water | 1.1774e-07 | 5.2071e-05 | -5.9104e-05 | 7.1467e-06 | -3.3560e-07 | 1.0797 |
| 1-hexylamine | -6.8038e-06 | -1.0323e-04 | 1.1422e-04 | -1.7932e-05 | 1.1438e-06 | 176.0360 |
| methanol | -1.7309e-06 | -3.1290e-06 | 3.4408e-06 | -2.0605e-06 | 6.4805e-06 | 4.1867 |
| ethene | -5.7397e-07 | 4.1682e-05 | -4.6175e-05 | 3.9049e-06 | 1.8240e-07 | 3.3631 |
| methane | 2.5309e-08 | -4.5112e-06 | 5.3238e-06 | -7.9313e-07 | 9.0110e-07 | 1.9833 |
| ethylamine | -7.6544e-06 | -4.9935e-05 | 5.5777e-05 | -1.3533e-05 | -9.1540e-07 | 18.0129 |
| ethanol | 1.1160e-05 | -4.8624e-05 | 5.6853e-05 | 2.8938e-06 | 7.5760e-06 | 14.6105 |
| hexane | 1.4781e-05 | 2.8967e-05 | -3.1717e-05 | 1.7421e-05 | -1.1650e-07 | 130.9111 |
| cyclopentane | 1.9080e-05 | 8.7994e-05 | -9.3962e-05 | 2.4965e-05 | 2.6800e-08 | 54.5885 |
| butane | 1.2717e-05 | -3.5731e-05 | 3.7339e-05 | 1.1050e-05 | 5.9610e-07 | 41.5918 |
| 1-pentanol | -5.7729e-05 | -9.9320e-05 | 1.2619e-04 | -8.4705e-05 | 5.0051e-06 | 102.8831 |
| 1-butanol | 2.6880e-05 | -9.4562e-05 | 1.0508e-04 | 1.6285e-05 | 9.2200e-07 | 74.2747 |
| cyclopropane | 1.1728e-05 | 4.1420e-05 | -4.5669e-05 | 1.5943e-05 | -9.2150e-07 | 11.6828 |
| hydrogen_peroxide | -2.2752e-07 | 2.0495e-05 | -2.3157e-05 | 2.4183e-06 | 2.3700e-07 | 2.3459 |
| ammonia | 3.7417e-07 | 2.9976e-05 | -3.3210e-05 | 3.6035e-06 | 9.6700e-08 | 1.6706 |
| cyclobutane | -1.0257e-05 | 1.7592e-04 | -2.0002e-04 | 1.3786e-05 | -1.0660e-05 | 22.9783 |
| pentane | -5.5838e-06 | -2.1463e-04 | 2.2873e-04 | -1.9761e-05 | 2.0828e-06 | 79.4211 |
| pyridine | 4.5146e-06 | -4.4350e-05 | 4.7036e-05 | 1.7380e-06 | 2.0220e-07 | 42.1721 |
| 1-propylamine | 1.1753e-05 | 6.6638e-05 | -7.4420e-05 | 1.9476e-05 | -3.9660e-07 | 51.0550 |
| 1-propene | 2.9772e-07 | 7.4772e-06 | -8.6289e-06 | 1.4183e-06 | -7.0850e-07 | 13.0076 |
| carbon_dioxide | 1.4071e-06 | -3.5163e-05 | 4.2448e-05 | -5.8988e-06 | -8.6200e-08 | 1.9148 |
| cyclohexane | -4.9942e-05 | 2.8359e-04 | -3.3186e-04 | -1.7792e-06 | -3.0933e-05 | 84.4095 |
| 1-pentene | -3.4967e-05 | -2.7088e-05 | 3.3169e-05 | -4.1124e-05 | -3.2000e-08 | 64.0157 |
| 1-pentylamine | -9.2611e-06 | -6.9553e-05 | 7.6401e-05 | -1.6219e-05 | -1.9798e-06 | 129.7257 |
| 1-butene | 1.1407e-05 | -1.2461e-04 | 1.5473e-04 | -1.8756e-05 | 1.4308e-06 | 38.8057 |
| ethane | 2.5769e-06 | -1.5856e-05 | 1.7729e-05 | 6.8581e-07 | -6.0010e-07 | 4.9346 |
| 1-hexene | 2.8032e-06 | -1.2182e-05 | 8.5036e-06 | 6.3801e-06 | 1.4920e-07 | 136.3154 |

- **BLYP**

| System | $\quad\Delta E_\text{tot}\quad$ | $\quad\Delta E_\text{1e}\quad$ | $\quad\Delta E_\text{2e}\quad$ | $\quad\Delta E_\text{xc}\quad$ | $\Delta E_\text{HOMO-LUMO gap}$ | Wall time (sec) | 
| :- | :-: | :-: | :-: | :-: | :-: | -: |
| methylamine | -3.4158e-06 | -5.4660e-06 | 5.3908e-06 | -3.3592e-06 | 2.2086e-06 | 8.1370 |
| carbon_monoxide | 7.6566e-07 | -1.1822e-04 | 1.3371e-04 | -1.4737e-05 | -2.4410e-07 | 2.8111 |
| benzene | 5.5895e-06 | 7.0633e-05 | -7.6796e-05 | 1.1663e-05 | 2.0710e-07 | 41.9015 |
| 1-hexanol | -3.7437e-05 | -1.7841e-05 | 2.7389e-05 | -4.7127e-05 | -7.6600e-07 | 156.7170 |
| oxirane | -9.4525e-07 | -7.0850e-06 | 7.8288e-06 | -1.7219e-06 | 5.9990e-07 | 12.8891 |
| pyrrole | -3.6845e-06 | 2.5182e-05 | -2.8764e-05 | -1.7326e-07 | 2.3940e-07 | 31.0708 |
| furan | -2.4140e-06 | 1.3709e-05 | -1.4085e-05 | -2.1086e-06 | -1.4460e-07 | 30.7187 |
| propane | -1.1635e-05 | -1.1520e-04 | 1.2569e-04 | -2.2161e-05 | -6.2900e-07 | 21.5090 |
| 1-propanol | 1.1320e-05 | 4.9515e-06 | -7.5846e-06 | 1.3894e-05 | -1.4827e-06 | 79.1883 |
| 1-butylamine | 4.9500e-07 | 2.9983e-05 | -3.2166e-05 | 2.5956e-06 | -2.3196e-06 | 77.7563 |
| water | 2.1403e-07 | -1.5819e-04 | 1.7886e-04 | -2.0468e-05 | 3.4300e-08 | 1.7797 |
| 1-hexylamine | -4.4126e-06 | -1.5205e-04 | 1.6478e-04 | -1.7279e-05 | 5.2580e-07 | 176.4170 |
| methanol | -1.7133e-06 | -6.6279e-06 | 8.2219e-06 | -3.3250e-06 | 4.2416e-06 | 6.9046 |
| ethene | -2.8329e-07 | 1.4407e-04 | -1.6367e-04 | 1.9302e-05 | -5.8570e-07 | 5.2840 |
| methane | 4.3539e-08 | -9.4489e-06 | 1.0604e-05 | -1.1180e-06 | 2.7670e-07 | 3.1808 |
| ethylamine | -7.4635e-06 | -7.0719e-05 | 7.9068e-05 | -1.5850e-05 | -4.7390e-07 | 21.2482 |
| ethanol | 1.3092e-05 | -6.8260e-05 | 7.7500e-05 | 3.8160e-06 | 6.7762e-06 | 18.5641 |
| hexane | 1.3382e-05 | 5.1310e-05 | -5.4904e-05 | 1.6867e-05 | -1.4446e-06 | 134.1697 |
| cyclopentane | 1.6050e-05 | 8.8102e-05 | -9.4682e-05 | 2.2547e-05 | 8.3010e-07 | 57.1636 |
| butane | 9.7959e-06 | 2.2874e-05 | -2.5640e-05 | 1.2503e-05 | 3.6290e-07 | 44.3323 |
| 1-pentanol | -4.8670e-05 | 1.3929e-05 | 2.9612e-06 | -6.5672e-05 | 3.2538e-06 | 111.4073 |
| 1-butanol | 2.7061e-05 | 6.8702e-05 | -6.9582e-05 | 2.7858e-05 | -9.7700e-07 | 84.2336 |
| cyclopropane | 9.2693e-06 | 4.0419e-05 | -4.4120e-05 | 1.2937e-05 | -3.2500e-07 | 14.8105 |
| hydrogen_peroxide | -4.3466e-07 | 7.9258e-05 | -8.9659e-05 | 9.9500e-06 | 8.2650e-07 | 4.2311 |
| ammonia | 4.0432e-07 | 1.5176e-05 | -1.6666e-05 | 1.8891e-06 | 3.6860e-07 | 2.7723 |
| cyclobutane | -1.1024e-05 | 1.3205e-04 | -1.5119e-04 | 8.0648e-06 | -3.4779e-06 | 28.4853 |
| pentane | -4.1732e-06 | -1.4490e-04 | 1.5469e-04 | -1.4052e-05 | 1.3412e-06 | 82.2843 |
| pyridine | 2.3877e-06 | -6.2271e-05 | 6.6459e-05 | -1.8908e-06 | 7.8800e-07 | 44.7745 |
| 1-propylamine | 1.0828e-05 | -9.5405e-05 | 9.8680e-05 | 7.4945e-06 | 1.3623e-06 | 57.1780 |
| 1-propene | -6.0543e-07 | 5.1962e-06 | -6.4260e-06 | 5.9320e-07 | -8.9520e-07 | 17.6234 |
| carbon_dioxide | 1.6545e-06 | -5.6646e-05 | 6.8060e-05 | -9.7816e-06 | 1.1220e-07 | 4.6161 |
| cyclohexane | -2.7313e-05 | 2.2525e-04 | -2.6977e-04 | 1.7092e-05 | -2.7817e-05 | 87.1533 |
| 1-pentene | -2.9387e-05 | 1.4313e-05 | -1.1801e-05 | -3.1975e-05 | -2.6790e-07 | 65.2221 |
| 1-pentylamine | -8.7238e-06 | 1.9610e-05 | -1.8562e-05 | -9.8809e-06 | -2.6528e-06 | 127.4520 |
| 1-butene | 1.4674e-05 | -1.1320e-04 | 1.4084e-04 | -1.3023e-05 | 2.1204e-06 | 42.7310 |
| ethane | 2.1856e-06 | -2.2522e-05 | 2.5187e-05 | -4.9816e-07 | -4.8170e-07 | 8.8190 |
| 1-hexene | 3.2542e-06 | 3.2270e-05 | -3.3780e-05 | 4.6620e-06 | -8.2770e-07 | 142.2336 |

- **HF**

| System | $\quad\Delta E_\text{tot}\quad$ | Wall time (sec) | 
| :- | :-: | -: |
| methanol | 3.1604e-10 | 3.003016 | 
| ethanol | 3.5021e-10 | 13.171072 | 
| 1-propanol | 6.8701e-10 | 35.912599 | 
| 1-butanol | 2.6955e-10 | 96.356437 | 
| 1-pentanol | -5.6957e-11 | 195.1718 | 
| 1-hexanol | 2.0833e-10 | 353.214047 | 
| methane | 4.5851e-11 | 0.945244 | 
| ethane | -9.1902e-11 | 5.400898 | 
| propane | -5.9757e-11 | 18.332309 | 
| butane | 2.2504e-10 | 5.672531 | 
| pentane | -5.0156e-10 | 124.540576 | 
| hexane | -4.2908e-10 | 252.264091 | 
| ethene | 2.9544e-10 | 3.192671 | 
| 1-propene | -1.7594e-10 | 11.814476 | 
| 1-butene | -5.4098e-10 | 39.080278 | 
| 1-pentene | -6.9781e-10 | 92.765139 | 
| 1-hexene | -9.8126e-10 | 187.272921 | 
| methylamine | 1.3522e-10 | 4.016299 | 
| ethylamine | 9.5213e-11 | 16.972886 | 
| 1-propylamine | 2.1083e-10 | 48.243879 | 
| 1-butylamine | 7.1347e-10 | 110.805998 | 
| 1-pentylamine | 2.4392e-10 | 219.96451 | 
| 1-hexylamine | -4.2195e-10 | 389.22106 | 
| cyclopropane | 8.4412e-11 | 12.852455 | 
| cyclobutane | 5.5962e-11 | 38.54993 | 
| cyclopentane | -1.4057e-10 | 91.709966 | 
| cyclohexane | 1.3378e-10 | 192.254076 | 
| water | 3.1559e-10 | 0.438969 | 
| ammonia | 3.7608e-10 | 0.632635 | 
| pyrrole | 5.8455e-10 | 40.915752 | 
| pyridine | 8.0149e-11 | 72.489318 | 
| oxirane | 8.2983e-10 | 7.919546 | 
| hydrogen_peroxide | -6.6203e-10 | 1.643334 | 
| furan | 1.1821e-09 | 35.161369 | 
| carbon_monoxide | 3.5479e-09 | 0.761816 | 
| carbon_dioxide | 5.7953e-09 | 2.485237 | 


## 贡献该项目
该项目欢迎任何形式的贡献，包括但不限于新功能添加，性能优化，bug修复，文档更新等。

## 未来规划
增加力的计算。

---

[go to top](#)
