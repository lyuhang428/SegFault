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

    include/dft.hpp
    include/hf.hpp
    include/mints.hpp
    include/beckegrid.hpp
    include/lebedev.hpp
    include/rsh.hpp
)

target_include_directories(segfault
PRIVATE
    /path/to/include/eigen3/
    /path/to/include/xtl/
    /path/to/include/xtensor/
    /path/to/include/xtensor-blas/
    /path/to/include/libint/
    /path/to/include/libxc/
    /path/to/include/gsl/
)

target_link_directories(segfault
PRIVATE
    /path/to/lib/libint/
    /path/to/lib/libxc/
    /path/to/lib/gsl/
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

