cmake_minimum_required(VERSION 3.23)
project(segfault LANGUAGES CXX VERSION 0.0.1)

set(VAR_TIMEIT "ON")
set(VAR_RELEASE "ON")

add_executable(segfault)

target_sources(segfault
PRIVATE
    main.cc
    src/dft.cc
    src/hf.cc
    src/mints.cc
    src/beckegrid.cc
    src/lebedev.cc

    include/dft.hpp
    include/hf.hpp
    include/mints.hpp
    include/beckegrid.hpp
    include/lebedev.hpp
    include/lblas.hpp
)

target_compile_features(segfault PRIVATE cxx_std_20)
target_compile_definitions(segfault 
PRIVATE 
    RELEASE=${VAR_RELEASE}
    TIMEIT=${VAR_TIMEIT}
)

message("MACRO definitions : ")
message("RELEASE is ${VAR_RELEASE}")
message("TIMEIT is ${VAR_TIMEIT}")

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU") # g++
    target_include_directories(segfault
    PRIVATE
        /path/to/intel/oneapi/mkl/version/include/
    )

    target_link_directories(segfault
    PRIVATE
        /path/to/intel/oneapi/mkl/version/lib/
    )

    target_compile_options(segfault PRIVATE -Wall -Wno-sign-compare -Wno-unused-variable -m64 -fopenmp) # how to link the right openmp lib?    
    target_link_options(segfault PRIVATE -m64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl)

elseif(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM") # icpx
    target_compile_options(segfault PRIVATE -Wall -Wno-unused-variable -qmkl -qopenmp)
    target_link_options(segfault PRIVATE -qmkl -qopenmp)
endif()


# path to xtensor header
# target_include_directories(segfault
# PRIVATE

# )


# lib to libxc and libint
# target_link_directories(segfault
# PRIVATE

# )


target_link_libraries(segfault
PRIVATE
    int2     # libint2
    xc       # libxc
)

install(TARGETS segfault)

