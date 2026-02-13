#!/usr/bin/bash
# just less typing...

installdir=install
builddir=build

if [[ -d "$installdir" ]] ; then echo "deleting ./install folder" ; rm -rf ${installdir} ; fi
if [[ -d "$builddir" ]] ; then echo "deleting ./build folder" ; rm -rf ${builddir} ; fi

export DATADIR=${PWD}/data
cmake -B build -DCMAKE_CXX_COMPILER=icpx -DCMAKE_INSTALL_PREFIX=./install
cmake --build build/ -j8
cmake --install build/
