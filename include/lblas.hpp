#pragma once 

#define Intel
// TODO: check if openblas and mkl use different arguments

#if defined GNU
#include "OpenBLAS/include/cblas.h"
#include "OpenBLAS/include/lapacke.h"
#elif defined Intel
#include "mkl_cblas.h"
#include "mkl_lapacke.h"
#endif

//>! compute trace
template<typename T>
inline T ltr(const T* data, const int n)
{
    T tr = 0;
    for (auto i=0; i < n; ++i) tr += *(data+n*i+i);
    return tr;
}


