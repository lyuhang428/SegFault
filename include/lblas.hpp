#ifndef INCLUDE_LBLAS_HPP
#define INCLUDE_LBLAS_HPP

#include "mkl_cblas.h"
#include "mkl_lapacke.h"

//>! Notes !<//
//>! 如果是行主序，lda使用的是矩阵A的列数！
//>! lda = leading dimension of a, in which dimension data in a is stored contiguously

//>! DSYEV(R/X) AX=XL will destroy original A

template<typename T>
inline T ltr(const T* data, const int n)
{
    T tr = 0;
    for (auto i=0; i < n; ++i) tr += *(data+n*i+i);
    return tr;
}

/*
void cblas_dgemm(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE TransA,
                 const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N,
                 const MKL_INT K, const double alpha, const double *A,
                 const MKL_INT lda, const double *B, const MKL_INT ldb,
                 const double beta, double *C, const MKL_INT ldc) NOTHROW;
//>! mxk @ kxn -> mxn
*/
// inline void ldgemm(const int nrowA, const int k, const int ncolB, const double* A, const double* B, double* C)
// {
//     // A nrowA x k
//     // B k x ncolB
//     // -> C nrowA x ncolB
//     cblas_dgemm(CblasRowMajor, CblasNoTrans, 
//                 CblasNoTrans, nrowA, ncolB, 
//                 k, 1., A, 
//                 k, B, ncolB, 
//                 0., C, ncolB); 
// }


/*
void cblas_dgemv(const CBLAS_LAYOUT Layout,
                  const CBLAS_TRANSPOSE TransA, const MKL_INT M, const MKL_INT N,
                  const double alpha, const double *A, const MKL_INT lda,
                  const double *X, const MKL_INT incX, const double beta,
                  double *Y, const MKL_INT incY) NOTHROW;
//>! mxn @ n -> m
*/
// inline void ldgemv(const int nrowA, const int ncolA, const double* A, const double* x, double* y)
// {
//     cblas_dgemv(CblasRowMajor, 
//                 CblasNoTrans, nrowA, ncolA, 
//                 1., A, ncolA, 
//                 x, 1, 0., 
//                 y, 1);
// }


/*
void cblas_dger(const CBLAS_LAYOUT Layout, const MKL_INT M, const MKL_INT N,
                const double alpha, const double *X, const MKL_INT incX,
                const double *Y, const MKL_INT incY, double *A, const MKL_INT lda) NOTHROW;
//>! x(m) outer y(n) -> A(mxn)
*/
// inline void ldger(const int m, const int n, const double* x, const double* y, double* A)
// {
//     cblas_dger(CblasRowMajor, m, n, 
//                 1., x, 1, 
//                 y, 1, A, n); // row-major, use ncol!
// }


/*
lapack_int LAPACKE_dsyev(int matrix_layout, char jobz, char uplo, lapack_int n,
                          double* a, lapack_int lda, double* w);
//>! A(nxn)
*/
// inline int ldsyev(const int n, double* vecs, double* vals)
// {
//     // A will be overwritten by its eigvecs. If jobz=eigvals, A will be desctroied
//     // 'V' compute vecs AND vals
//     // 'U': upper+diag of A; 'L': lower+diag of A
//     int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', n, 
//                              vecs, n, vals);
//    return info; 
// }


/*
lapack_int LAPACKE_dsyevr( int matrix_layout, char jobz, char range, char uplo,
                           lapack_int n, double* a, lapack_int lda, double vl,
                           double vu, lapack_int il, lapack_int iu,
                           double abstol, lapack_int* m, double* w, double* z,
                           lapack_int ldz, lapack_int* isuppz );
//>! A(nxn) -> vecs(nxm)
*/
// inline int ldsyevr(const int n, const int nvals, double* A, double* vecs, double* vals)
// {
//     // A will be destroied anyway no matter what jobz
//     // w(n), first nvals are requested eigvals, others not affected
//     // z(nxn), first nvals col are requested eigvecs, others not affected
//     double sfmin = LAPACKE_dlamch('S'); // safe minium
//     int m; 
//     std::vector<int> isuppz(2*nvals, 0);
    
//     int info = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, 'V', 'I', 'U', 
//                               n, A, n, 0., 
//                               0., 1, nvals, 
//                               sfmin, &m, vals, vecs, 
//                               n, isuppz.data());
                              
//     return info;
// }

//>! similar to dsyevr, to be interfaced
// inline int ldsyevx();


/*
lapack_int LAPACKE_dsytrf( int matrix_layout, char uplo, lapack_int n, double* a,
                           lapack_int lda, lapack_int* ipiv );

lapack_int LAPACKE_dsytri( int matrix_layout, char uplo, lapack_int n, double* a,
                           lapack_int lda, const lapack_int* ipiv );
//>! dsytrf(A) -> dsytri(A=U**T*D*U|L*D*L**T) -> A^-1; A is real symmetric
*/
// inline int linv(const int n, double* A)
// {
//     // A is symmetric, only upper triangle part of A is overwritten by its inverse
//     std::vector<int> ipiv(n, 0);
//     int info = LAPACKE_dsytrf(LAPACK_ROW_MAJOR, 'U', n, A, 
//         n, ipiv.data()); // A -> U**T*D*U|L*D*L**T
        
//     assert(info == 0);
//     info = LAPACKE_dsytri(LAPACK_ROW_MAJOR, 'U', n, A, 
//                             n, ipiv.data()); // A -> A^-1
            
//     return info;
// }


/*
lapack_int LAPACKE_dgesv( int matrix_layout, lapack_int n, lapack_int nrhs,
                          double* a, lapack_int lda, lapack_int* ipiv,
                          double* b, lapack_int ldb );
//>! Ax = b ; A is real general
*/
// inline int ldgesv(const int n, double* A, double* b)
// {
//     std::vector<int> ipiv(n, 0);
//     int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, 1, 
//                               A, n, ipiv.data(), 
//                               b, 1);
//     return info;
// }


/*
lapack_int LAPACKE_dsysv( int matrix_layout, char uplo, lapack_int n,
                          lapack_int nrhs, double* a, lapack_int lda,
                          lapack_int* ipiv, double* b, lapack_int ldb );
//>! Ax = b ; A is real symmetric
*/
// inline int ldsysv(const int n, double* A, double* b)
// {
//     std::vector<int> ipiv(n, 0);
//     int info = LAPACKE_dsysv(LAPACK_ROW_MAJOR, 'U', n, 
//                             1, A, n, 
//                             ipiv.data(), b, 1);
//     return info;
// }


#endif
