#!/usr/bin/bash
# just less typing...

installdir=install
builddir=build

#if [[ -d "$installdir" ]] ; then echo "deleting ./install folder" ; rm -rf ${installdir} ; fi
#if [[ -d "$builddir" ]] ; then echo "deleting ./build folder" ; rm -rf ${builddir} ; fi

is_intel=$(whereis icpx | cut -d ':' -f 2)
if [[ -z $is_intel ]] ; then echo 'source intel setvars' ; source /path/to/intel/oneapi/setvars.sh ; fi

cmake -B build -DCMAKE_CXX_COMPILER=icpx -DCMAKE_INSTALL_PREFIX=./install
cmake --build build/ -j4
cmake --install build/
