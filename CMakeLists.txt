cmake_minimum_required(VERSION 3.23)
project(segfault LANGUAGES CXX VERSION 0.0.0)

add_executable(segfault) # add_library() if want to build lib

target_sources(segfault
PRIVATE
    main.cc
    src/dft.cc
    src/hf.cc
    src/mints.cc
    src/beckegrid.cc
    src/lebedev.cc
    src/rsh.cc

PUBLIC
    include/dft.hpp
    include/hf.hpp
    include/mints.hpp
    include/beckegrid.hpp
    include/lebedev.hpp
    include/rsh.hpp
)

target_compile_features(segfault PUBLIC cxx_std_20)
target_compile_options(segfault PUBLIC -Wall -qmkl -qopenmp)
# target_compile_definitions() # 等价于 -D 宏定义
target_link_options(segfault PUBLIC -qmkl -qopenmp) # mkl omp needs link option too

# on PATH, so no need to add_library INTERFACE
target_link_libraries(segfault
PUBLIC
    int2     # libint2
    xc       # libxc
    gsl      # libgsl
    gslcblas # libgslcblas
)

install(TARGETS segfault)

# 通过添加目录可以链接非PATH上的库？
# add_library(VendorLib INTERFACE) # INTERFACE means no build, 3rd-party library
# target_include_directories(VendorLib INTERFACE include)
# target_link_directories(VendorLib INTERFACE lib)
# target_link_libraries(VendorLib INTERFACE Vendor)

# 上一级CMake
# target_link_libraries() # link libVendorLib
# 再上一级
# add_subdirectory()
