#include "RBSim/lcp/matrix.h"

#include <cstring>

// misc defines

#define dMACRO_MAX(a, b) ((a) > (b) ? (a) : (b))
#define dMACRO_MIN(a, b) ((a) < (b) ? (a) : (b))

void dxMultiply0(double *A, const double *B, const double *C, unsigned p, unsigned q, unsigned r)
{
    assert (A && B && C && p>0 && q>0 && r>0);
    const unsigned qskip = dPAD(q);
    const unsigned rskip = dPAD(r);
    double *aa = A;
    const double *bb = B;
    for (unsigned i = p; i != 0; aa+=rskip, bb+=qskip, --i) {
        double *a = aa;
        const double *cc = C, *ccend = C + r;
        for (; cc != ccend; ++a, ++cc) {
            double sum = 0.0;
            const double *c = cc;
            const double *b = bb, *bend = bb + q;
            for (; b != bend; c+=rskip, ++b) {
                sum += (*b) * (*c);
            }
            (*a) = sum;
        }
    }
}


/*extern */
void dxMultiply1(double *A, const double *B, const double *C, unsigned p, unsigned q, unsigned r)
{
    assert (A && B && C && p>0 && q>0 && r>0);
    const unsigned pskip = dPAD(p);
    const unsigned rskip = dPAD(r);
    double *aa = A;
    const double *bb = B, *bbend = B + p;
    for (; bb != bbend; aa += rskip, ++bb) {
        double *a = aa;
        const double *cc = C, *ccend = C + r;
        for (; cc != ccend; ++a, ++cc) {
            double sum = 0.0;
            const double *b = bb, *c = cc;
            for (unsigned k = q; k != 0; b += pskip, c += rskip, --k) {
                sum += (*b) * (*c);
            }
            (*a) = sum;
        }
    }
}


/*extern */
void dxMultiply2(double *A, const double *B, const double *C, unsigned p, unsigned q, unsigned r)
{
    assert (A && B && C && p>0 && q>0 && r>0);
    const unsigned rskip = dPAD(r);
    const unsigned qskip = dPAD(q);
    double *aa = A;
    const double *bb = B;
    for (unsigned i = p; i != 0; aa += rskip, bb += qskip, --i) {
        double *a = aa, *aend = aa + r;
        const double *cc = C;
        for (; a != aend; cc+=qskip, ++a) {
            double sum = 0.0;
            const double *b = bb, *c = cc, *cend = cc + q;
            for (; c != cend; ++b, ++c) {
                sum += (*b) * (*c);
            }
            (*a) = sum;
        }
    }
}


/*extern */
int dxFactorCholesky(double *A, unsigned n, void *tmpBuf/*[n]*/)
{
    assert (n > 0 && A);
    bool failure = false;

    double *alloctedBuf = NULL;
    size_t allocatedSize;

    const unsigned nskip = dPAD (n);

    double *recip = (double *)tmpBuf;
    if (tmpBuf == NULL) {
        allocatedSize = n * sizeof(double);
        alloctedBuf = allocatedSize > STACK_ALLOC_MAX ? (double *)malloc(allocatedSize) : NULL;
        recip = alloctedBuf != NULL ? alloctedBuf : (double*)ALLOCA(allocatedSize);
    }

    double *aa = A;
    for (unsigned i = 0; i < n; aa += nskip, ++i) {
        double *cc = aa;
        {
            const double *bb = A;
            for (unsigned j = 0; j < i; bb += nskip, ++cc, ++j) {
                double sum = *cc;
                const double *a = aa, *b = bb, *bend = bb + j;
                for (; b != bend; ++a, ++b) {
                    sum -= (*a) * (*b);
                }
                *cc = sum * recip[j];
            }
        }
        {
            double sum = *cc;
            double *a = aa, *aend = aa + i;
            for (; a != aend; ++a) {
                sum -= (*a)*(*a);
            }
            if (sum <= 0.0) {
                failure = true;
                break;
            }
            double sumsqrt = sqrt(sum);
            *cc = sumsqrt;
            recip[i] = 1.0 / (sumsqrt);
        }
    }

    if (alloctedBuf != NULL) {
        free(alloctedBuf);
    }

    return failure ? 0 : 1;
}

/*extern */
void dxSolveCholesky(const double *L, double *b, unsigned n, void *tmpBuf/*[n]*/)
{
    assert (n > 0 && L && b);

    double *alloctedBuf = NULL;
    size_t allocatedSize;

    const unsigned nskip = dPAD (n);

    double *y = (double *)tmpBuf;
    if (tmpBuf == NULL) {
        allocatedSize = n * sizeof(double);
        alloctedBuf = allocatedSize > STACK_ALLOC_MAX ? (double *)malloc(allocatedSize) : NULL;
        y = alloctedBuf != NULL ? alloctedBuf : (double*)ALLOCA(allocatedSize);
    }

    {
        const double *ll = L;
        for (unsigned i = 0; i < n; ll += nskip, ++i) {
            double sum = 0.0;
            for (unsigned k = 0; k < i; ++k) {
                sum += ll[k] * y[k];
            }
            assert(ll[i] != 0.0);
            y[i] = (b[i] - sum) / ll[i];
        }
    }
    {
        const double *ll = L + (n - 1) * (nskip + 1);
        for (unsigned i = n; i > 0; ll -= nskip + 1) {
            --i;
            double sum = 0.0;
            const double *l = ll + nskip;
            for (unsigned k = i + 1; k < n; l += nskip, ++k) {
                sum += (*l) * b[k];
            }
            assert(*ll != 0.0);
            b[i] = (y[i] - sum) / (*ll);
        }
    }

    if (alloctedBuf != NULL) {
        free(alloctedBuf);
    }
}


/*extern */
int dxInvertPDMatrix(const double *A, double *Ainv, unsigned n, void *tmpBuf/*[nskip*(n+2)]*/)
{
    assert (n > 0 && A && Ainv);
    bool success = false;

    double *alloctedBuf = NULL;
    size_t allocatedSize;

    size_t choleskyFactorSize = dxEstimateFactorCholeskyTmpbufSize(n);
    size_t choleskySolveSize = dxEstimateSolveCholeskyTmpbufSize(n);
    size_t choleskyMaxSize = dMACRO_MAX(choleskyFactorSize, choleskySolveSize);
    assert(choleskyMaxSize % sizeof(double) == 0);

    const unsigned nskip = dPAD (n);
    const size_t nskip_mul_n = (size_t)nskip * n;

    double *tmp = (double *)tmpBuf;
    if (tmpBuf == NULL) {
        allocatedSize = choleskyMaxSize + (nskip + nskip_mul_n) * sizeof(double);
        alloctedBuf = allocatedSize > STACK_ALLOC_MAX ? (double *)malloc(allocatedSize) : NULL;
        tmp = alloctedBuf != NULL ? alloctedBuf : (double*)ALLOCA(allocatedSize);
    }

    double *X = (double *)((char *)tmp + choleskyMaxSize);
    double *L = X + nskip;
    memcpy (L, A, nskip_mul_n * sizeof(double));
    if (dxFactorCholesky(L, n, tmp)) {
        dSetZero (Ainv, nskip_mul_n);	// make sure all padding elements set to 0
        double *aa = Ainv, *xi = X, *xiend = X + n;
        for (; xi != xiend; ++aa, ++xi) {
            dSetZero(X, n);
            *xi = 1.0;
            dxSolveCholesky(L, X, n, tmp);
            double *a = aa;
            const double *x = X, *xend = X + n;
            for (; x != xend; a += nskip, ++x) {
                *a = *x;
            }
        }
        success = true;
    }

    if (alloctedBuf != NULL) {
        free(alloctedBuf);
    }

    return success ? 1 : 0;
}


/*extern */
int dxIsPositiveDefinite(const double *A, unsigned n, void *tmpBuf/*[nskip*(n+1)]*/)
{
    assert (n > 0 && A);

    double *alloctedBuf = NULL;
    size_t allocatedSize;

    size_t choleskyFactorSize = dxEstimateFactorCholeskyTmpbufSize(n);
    assert(choleskyFactorSize % sizeof(double) == 0);

    const unsigned nskip = dPAD (n);
    const size_t nskip_mul_n = (size_t)nskip * n;

    double *tmp = (double *)tmpBuf;
    if (tmpBuf == NULL) {
        allocatedSize = choleskyFactorSize + nskip_mul_n * sizeof(double);
        alloctedBuf = allocatedSize > STACK_ALLOC_MAX ? (double *)malloc(allocatedSize) : NULL;
        tmp = alloctedBuf != NULL ? alloctedBuf : (double*)ALLOCA(allocatedSize);
    }

    double *Acopy = (double *)((char *)tmp + choleskyFactorSize);
    memcpy(Acopy, A, nskip_mul_n * sizeof(double));
    int factorResult = dxFactorCholesky (Acopy, n, tmp);

    if (alloctedBuf != NULL) {
        free(alloctedBuf);
    }

    return factorResult;
}


/*extern */
void dxLDLTAddTL(double *L, double *d, const double *a, unsigned n, unsigned nskip, void *tmpBuf/*[2*nskip]*/)
{
    assert(L && d && a && n > 0 && nskip >= n);

    if (n < 2) return;

    double *alloctedBuf = NULL;
    size_t allocatedSize;

    double *W1 = (double *)tmpBuf;
    if (tmpBuf == NULL) {
        allocatedSize = nskip * (2 * sizeof(double));
        alloctedBuf = allocatedSize > STACK_ALLOC_MAX ? (double *)malloc(allocatedSize) : NULL;
        W1 = alloctedBuf != NULL ? alloctedBuf : (double*)ALLOCA(allocatedSize);
    }

    double *W2 = W1 + nskip;

    W1[0] = 0.0;
    W2[0] = 0.0;
    for (unsigned j = 1; j < n; ++j) {
        W1[j] = W2[j] = (double) (a[j] * M_SQRT1_2);
    }
    double W11 = (double) ((0.5*a[0]+1)*M_SQRT1_2);
    double W21 = (double) ((0.5*a[0]-1)*M_SQRT1_2);

    double alpha1 = 1.0;
    double alpha2 = 1.0;

    {
        double dee = d[0];
        double alphanew = alpha1 + (W11*W11)*dee;
        assert(alphanew != double(0.0));
        dee /= alphanew;
        double gamma1 = W11 * dee;
        dee *= alpha1;
        alpha1 = alphanew;
        alphanew = alpha2 - (W21*W21)*dee;
        dee /= alphanew;
        //double gamma2 = W21 * dee;
        alpha2 = alphanew;
        double k1 = 1.0 - W21*gamma1;
        double k2 = W21*gamma1*W11 - W21;
        double *ll = L + nskip;
        for (unsigned p = 1; p < n; ll += nskip, ++p) {
            double Wp = W1[p];
            double ell = *ll;
            W1[p] =    Wp - W11*ell;
            W2[p] = k1*Wp +  k2*ell;
        }
    }

    double *ll = L + (nskip + 1);
    for (unsigned j = 1; j < n; ll += nskip + 1, ++j) {
        double k1 = W1[j];
        double k2 = W2[j];

        double dee = d[j];
        double alphanew = alpha1 + (k1*k1)*dee;
        assert(alphanew != double(0.0));
        dee /= alphanew;
        double gamma1 = k1 * dee;
        dee *= alpha1;
        alpha1 = alphanew;
        alphanew = alpha2 - (k2*k2)*dee;
        dee /= alphanew;
        double gamma2 = k2 * dee;
        dee *= alpha2;
        d[j] = dee;
        alpha2 = alphanew;

        double *l = ll + nskip;
        for (unsigned p = j + 1; p < n; l += nskip, ++p) {
            double ell = *l;
            double Wp = W1[p] - k1 * ell;
            ell += gamma1 * Wp;
            W1[p] = Wp;
            Wp = W2[p] - k2 * ell;
            ell -= gamma2 * Wp;
            W2[p] = Wp;
            *l = ell;
        }
    }

    if (alloctedBuf != NULL) {
        free(alloctedBuf);
    }
}


// macros for dLDLTRemove() for accessing A - either access the matrix
// directly or access it via row pointers. we are only supposed to reference
// the lower triangle of A (it is symmetric), but indexes i and j come from
// permutation vectors so they are not predictable. so do a test on the
// indexes - this should not slow things down too much, as we don't do this
// in an inner loop.

#define _GETA(i,j) (A[i][j])
//#define _GETA(i,j) (A[(i)*nskip+(j)])
#define GETA(i,j) ((i > j) ? _GETA(i,j) : _GETA(j,i))


/*extern */
void dxLDLTRemove(double **A, const unsigned *p, double *L, double *d,
    unsigned n1, unsigned n2, unsigned r, unsigned nskip, void *tmpBuf/*n2 + 2*nskip*/)
{
    assert(A && p && L && d && n1 > 0 && n2 > 0 /*&& r >= 0 */&& r < n2 &&
        n1 >= n2 && nskip >= n1);
#ifndef dNODEBUG
    for (unsigned i = 0; i < n2; ++i) assert(p[i] >= 0 && p[i] < n1);
#endif

    if (r == n2 - 1) {
        return;		// deleting the last row/col is easy
    }

    double *alloctedBuf = NULL;
    size_t allocatedSize;

    size_t LDLTAddTLSize = dxEstimateLDLTAddTLTmpbufSize(nskip);
    assert(LDLTAddTLSize % sizeof(double) == 0);

    double *tmp = (double *)tmpBuf;
    if (tmpBuf == NULL) {
        allocatedSize = LDLTAddTLSize + n2 * sizeof(double);
        alloctedBuf = allocatedSize > STACK_ALLOC_MAX ? (double *)malloc(allocatedSize) : NULL;
        tmp = alloctedBuf != NULL ? alloctedBuf : (double*)ALLOCA(allocatedSize);
    }
    
    if (r == 0) {
        double *a = (double *)((char *)tmp + LDLTAddTLSize);
        const unsigned p_0 = p[0];
        for (unsigned i = 0; i < n2; ++i) {
            a[i] = -GETA(p[i],p_0);
        }
        a[0] += 1.0;
        dxLDLTAddTL (L, d, a, n2, nskip, tmp);
    }
    else {
        double *t = (double *)((char *)tmp + LDLTAddTLSize);
        {
            double *Lcurr = L + r*nskip;
            for (unsigned i = 0; i < r; ++Lcurr, ++i) {
                assert(d[i] != double(0.0));
                t[i] = *Lcurr / d[i];
            }
        }
        double *a = t + r;
        {
            double *Lcurr = L + r * nskip;
            const unsigned *pp_r = p + r, p_r = *pp_r;
            const unsigned n2_minus_r = n2 - r;
            for (unsigned i = 0; i < n2_minus_r; Lcurr += nskip, ++i) {
                a[i] = dDot(Lcurr, t, r) - GETA(pp_r[i], p_r);
            }
        }
        a[0] += 1.0;
        dxLDLTAddTL (L + (size_t)(nskip + 1) * r, d + r, a, n2 - r, nskip, tmp);
    }

    // snip out row/column r from L and d
    dxRemoveRowCol (L, n2, nskip, r);
    if (r < (n2 - 1)) memmove (d + r, d + r + 1, (n2 - r - 1) * sizeof(double));

    if (alloctedBuf != NULL) {
        free(alloctedBuf);
    }
}


/*extern */
void dxRemoveRowCol(double *A, unsigned n, unsigned nskip, unsigned r)
{
    assert(A && n > 0 && nskip >= n && r >= 0 && r < n);
    if (r >= n - 1) return;
    if (r > 0) {
        {
            const size_t move_size = (n - r - 1) * sizeof(double);
            double *Adst = A + r;
            for (unsigned i = 0; i < r; Adst += nskip, ++i) {
                double *Asrc = Adst + 1;
                memmove (Adst, Asrc, move_size);
            }
        }
        {
            const size_t cpy_size = r * sizeof(double);
            double *Adst = A + (size_t)nskip * r;
            unsigned n1 = n - 1;
            for (unsigned i = r; i < n1; ++i) {
                double *Asrc = Adst + nskip;
                memcpy (Adst, Asrc, cpy_size);
                Adst = Asrc;
            }
        }
    }
    {
        const size_t cpy_size = (n - r - 1) * sizeof(double);
        double *Adst = A + (size_t)(nskip + 1) * r;
        unsigned n1 = n - 1;
        for (unsigned i = r; i < n1; ++i) {
            double *Asrc = Adst + (nskip + 1);
            memcpy (Adst, Asrc, cpy_size);
            Adst = Asrc - 1;
        }
    }
}


#undef dSetZero
#undef dSetValue
#undef dMultiply0
#undef dMultiply1
#undef dMultiply2
#undef dFactorCholesky
#undef dSolveCholesky
#undef dInvertPDMatrix
#undef dIsPositiveDefinite
#undef dLDLTAddTL
#undef dLDLTRemove
#undef dRemoveRowCol

/*extern ODE_API */
void dSetZero(double *a, int n)
{
    dxSetZero(a, n);
}

/*extern ODE_API */
void dSetValue(double *a, int n, double value)
{
    dxSetValue(a, n, value);
}

// double dDot (const double *a, const double *b, int n);

/*extern ODE_API */
void dMultiply0(double *A, const double *B, const double *C, int p,int q,int r)
{
    dxMultiply0(A, B, C, p, q, r);
}

/*extern ODE_API */
void dMultiply1(double *A, const double *B, const double *C, int p,int q,int r)
{
    dxMultiply1(A, B, C, p, q, r);
}

/*extern ODE_API */
void dMultiply2(double *A, const double *B, const double *C, int p,int q,int r)
{
    dxMultiply2(A, B, C, p, q, r);
}

/*extern ODE_API */
int dFactorCholesky(double *A, int n)
{
    return dxFactorCholesky(A, n, NULL);
}

/*extern ODE_API */
void dSolveCholesky(const double *L, double *b, int n)
{
    dxSolveCholesky(L, b, n, NULL);
}

/*extern ODE_API */
int dInvertPDMatrix (const double *A, double *Ainv, int n)
{
    return dxInvertPDMatrix(A, Ainv, n, NULL);
}

/*extern ODE_API */
int dIsPositiveDefinite(const double *A, int n)
{
    return dxIsPositiveDefinite(A, n, NULL);
}


/*extern ODE_API */
void dLDLTAddTL(double *L, double *d, const double *a, int n, int nskip)
{
    dxLDLTAddTL(L, d, a, n, nskip, NULL);
}

/*extern ODE_API */
void dLDLTRemove(double **A, const int *p, double *L, double *d, int n1, int n2, int r, int nskip)
{
    dxLDLTRemove(A, (const unsigned *)p, L, d, n1, n2, r, nskip, NULL);
    assert(sizeof(unsigned) == sizeof(*p));
}

/*extern ODE_API */
void dRemoveRowCol(double *A, int n, int nskip, int r)
{
    dxRemoveRowCol(A, n, nskip, r);
}

