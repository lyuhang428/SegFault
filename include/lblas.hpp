#ifndef INCLUDE_LBLAS_HPP
#define INCLUDE_LBLAS_HPP

#include "mkl_cblas.h"
#include "mkl_lapacke.h"

//>! this header includes lapacke and cblas
//>! TODO: add openblas support
//>! compute trace
template<typename T>
inline T ltr(const T* data, const int n)
{
    T tr = 0;
    for (auto i=0; i < n; ++i) tr += *(data+n*i+i);
    return tr;
}



#endif
