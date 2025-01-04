
/*
 * optimized and unoptimized vector and matrix functions 
 * (inlined private versions)
 */

#pragma once

#include <cmath>
#include <stddef.h>
#include <assert.h>

#include <RBSim/lcp/memory.h>

template <unsigned a_stride, typename element_type>
__inline
void dxtSetZero (element_type *a, size_t n)
{
    element_type *const aend = a + n * a_stride;
    for (element_type *acurr = a; acurr != aend; acurr += a_stride) {
        *acurr = (element_type)0;
    }
}

template <typename element_type>
__inline
void dxSetZero (element_type *a, size_t n) {
    dxtSetZero<1>(a, n);
}

template <typename element_type>
__inline
void dxSetValue (element_type *a, size_t n, element_type value) {
    element_type *const aend = a + n;
    for (element_type *acurr = a; acurr != aend; ++acurr) {
        *acurr = value;
    }
}

static __inline
double dxCalculateModuloMaximum(const double *valueStorage, size_t valueCount)
{
    assert(valueCount != 0);

    double moduleMaximum = fabs(*valueStorage);

    const double *const storageEnd = valueStorage + valueCount;
    for (const double *currentValue = valueStorage + 1; currentValue != storageEnd; ++currentValue) {
        moduleMaximum = fmax(moduleMaximum, fabs(*currentValue));
    }

    return moduleMaximum;
}


double dxDot (const double *a, const double *b, unsigned n);
void dxMultiply0 (double *A, const double *B, const double *C, unsigned p, unsigned q, unsigned r);
void dxMultiply1 (double *A, const double *B, const double *C, unsigned p, unsigned q, unsigned r);
void dxMultiply2 (double *A, const double *B, const double *C, unsigned p, unsigned q, unsigned r);
int dxFactorCholesky (double *A, unsigned n, void *tmpbuf);
void dxSolveCholesky (const double *L, double *b, unsigned n, void *tmpbuf);
int dxInvertPDMatrix (const double *A, double *Ainv, unsigned n, void *tmpbuf);
int dxIsPositiveDefinite (const double *A, unsigned n, void *tmpbuf);
void dxLDLTAddTL (double *L, double *d, const double *a, unsigned n, unsigned nskip, void *tmpbuf);
void dxLDLTRemove (double **A, const unsigned *p, double *L, double *d, unsigned n1, unsigned n2, unsigned r, unsigned nskip, void *tmpbuf);
void dxRemoveRowCol (double *A, unsigned n, unsigned nskip, unsigned r);

static __inline size_t dxEstimateFactorCholeskyTmpbufSize(unsigned n)
{
    return dPAD(n) * sizeof(double);
}

static __inline size_t dxEstimateSolveCholeskyTmpbufSize(unsigned n)
{
    return dPAD(n) * sizeof(double);
}

static __inline size_t dxEstimateInvertPDMatrixTmpbufSize(unsigned n)
{
    size_t FactorCholesky_size = dxEstimateFactorCholeskyTmpbufSize(n);
    size_t SolveCholesky_size = dxEstimateSolveCholeskyTmpbufSize(n);
    size_t MaxCholesky_size = FactorCholesky_size > SolveCholesky_size ? FactorCholesky_size : SolveCholesky_size;
    return (size_t)dPAD(n) * (n + 1) * sizeof(double) + MaxCholesky_size;
}

static __inline size_t dxEstimateIsPositiveDefiniteTmpbufSize(unsigned n)
{
    return (size_t)dPAD(n) * n * sizeof(double) + dxEstimateFactorCholeskyTmpbufSize(n);
}

static __inline size_t dxEstimateLDLTAddTLTmpbufSize(unsigned nskip)
{
    return nskip * (2 * sizeof(double));
}

static __inline size_t dxEstimateLDLTRemoveTmpbufSize(unsigned n2, unsigned nskip)
{
    return n2 * sizeof(double) + dxEstimateLDLTAddTLTmpbufSize(nskip);
}

/* For internal use */
#define dSetZero(a, n) dxSetZero(a, n)
#define dSetValue(a, n, value) dxSetValue(a, n, value)
#define dDot(a, b, n) dxDot(a, b, n)
#define dMultiply0(A, B, C, p, q, r) dxMultiply0(A, B, C, p, q, r)
#define dMultiply1(A, B, C, p, q, r) dxMultiply1(A, B, C, p, q, r)
#define dMultiply2(A, B, C, p, q, r) dxMultiply2(A, B, C, p, q, r)
#define dFactorCholesky(A, n, tmpbuf) dxFactorCholesky(A, n, tmpbuf)
#define dSolveCholesky(L, b, n, tmpbuf) dxSolveCholesky(L, b, n, tmpbuf)
#define dInvertPDMatrix(A, Ainv, n, tmpbuf) dxInvertPDMatrix(A, Ainv, n, tmpbuf)
#define dIsPositiveDefinite(A, n, tmpbuf) dxIsPositiveDefinite(A, n, tmpbuf)
#define dLDLTAddTL(L, d, a, n, nskip, tmpbuf) dxLDLTAddTL(L, d, a, n, nskip, tmpbuf)
#define dLDLTRemove(A, p, L, d, n1, n2, r, nskip, tmpbuf) dxLDLTRemove(A, p, L, d, n1, n2, r, nskip, tmpbuf)
#define dRemoveRowCol(A, n, nskip, r) dxRemoveRowCol(A, n, nskip, r)


#define dEstimateFactorCholeskyTmpbufSize(n) dxEstimateFactorCholeskyTmpbufSize(n)
#define dEstimateSolveCholeskyTmpbufSize(n) dxEstimateSolveCholeskyTmpbufSize(n)
#define dEstimateInvertPDMatrixTmpbufSize(n) dxEstimateInvertPDMatrixTmpbufSize(n)
#define dEstimateIsPositiveDefiniteTmpbufSize(n) dxEstimateIsPositiveDefiniteTmpbufSize(n)
#define dEstimateLDLTAddTLTmpbufSize(nskip) dxEstimateLDLTAddTLTmpbufSize(nskip)
#define dEstimateLDLTRemoveTmpbufSize(n2, nskip) dxEstimateLDLTRemoveTmpbufSize(n2, nskip)

