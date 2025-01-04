/*

THE ALGORITHM
-------------

solve A*x = b+w, with x and w subject to certain LCP conditions.
each x(i),w(i) must lie on one of the three line segments in the following
diagram. each line segment corresponds to one index set :

     w(i)
     /|\      |           :
      |       |           :
      |       |i in N     :
  w>0 |       |state[i]=0 :
      |       |           :
      |       |           :  i in C
  w=0 +       +-----------------------+
      |                   :           |
      |                   :           |
  w<0 |                   :           |i in N
      |                   :           |state[i]=1
      |                   :           |
      |                   :           |
      +-------|-----------|-----------|----------> x(i)
             lo           0           hi

the Dantzig algorithm proceeds as follows:
  for i=1:n
    * if (x(i),w(i)) is not on the line, push x(i) and w(i) positive or
      negative towards the line. as this is done, the other (x(j),w(j))
      for j<i are constrained to be on the line. if any (x,w) reaches the
      end of a line segment then it is switched between index sets.
    * i is added to the appropriate index set depending on what line segment
      it hits.

we restrict lo(i) <= 0 and hi(i) >= 0. this makes the algorithm a bit
simpler, because the starting point for x(i),w(i) is always on the dotted
line x=0 and x will only ever increase in one direction, so it can only hit
two out of the three line segments.


NOTES
-----

this is an implementation of "lcp_dantzig2_ldlt.m" and "lcp_dantzig_lohi.m".
the implementation is split into an LCP problem object (dLCP) and an LCP
driver function. most optimization occurs in the dLCP object.

a naive implementation of the algorithm requires either a lot of data motion
or a lot of permutation-array lookup, because we are constantly re-ordering
rows and columns. to avoid this and make a more optimized algorithm, a
non-trivial data structure is used to represent the matrix A (this is
implemented in the fast version of the dLCP object).

during execution of this algorithm, some indexes in A are clamped (set C),
some are non-clamped (set N), and some are "don't care" (where x=0).
A,x,b,w (and other problem vectors) are permuted such that the clamped
indexes are first, the unclamped indexes are next, and the don't-care
indexes are last. this permutation is recorded in the array `p'.
initially p = 0..n-1, and as the rows and columns of A,x,b,w are swapped,
the corresponding elements of p are swapped.

because the C and N elements are grouped together in the rows of A, we can do
lots of work with a fast dot product function. if A,x,etc were not permuted
and we only had a permutation array, then those dot products would be much
slower as we would have a permutation array lookup in some inner loops.

A is accessed through an array of row pointers, so that element (i,j) of the
permuted matrix is A[i][j]. this makes row swapping fast. for column swapping
we still have to actually move the data.

during execution of this algorithm we maintain an L*D*L' factorization of
the clamped submatrix of A (call it `AC') which is the top left nC*nC
submatrix of A. there are two ways we could arrange the rows/columns in AC.

(1) AC is always permuted such that L*D*L' = AC. this causes a problem
when a row/column is removed from C, because then all the rows/columns of A
between the deleted index and the end of C need to be rotated downward.
this results in a lot of data motion and slows things down.
(2) L*D*L' is actually a factorization of a *permutation* of AC (which is
itself a permutation of the underlying A). this is what we do - the
permutation is recorded in the vector C. call this permutation A[C,C].
when a row/column is removed from C, all we have to do is swap two
rows/columns and manipulate C.

*/

#include <cfloat>
#include <cmath>

#include "RBSim/lcp/lcp.h"

#include <algorithm>

#include "RBSim/lcp/fastdot_impl.h"
#include "RBSim/lcp/fastldltfactor_impl.h"
#include "RBSim/lcp/fastldltsolve_impl.h"
#include "RBSim/lcp/matrix.h"

#include "RBSim/lcp/memory.h"

//***************************************************************************
// code generation parameters

#define dLCP_FAST		// use fast dLCP object

#define NUB_OPTIMIZATIONS // use NUB optimizations

// option 1 : matrix row pointers (less data copying)
#define ROWPTRS
#define ATYPE double **
#define AROW(i) (m_A[i])

// option 2 : no matrix row pointers (slightly faster inner loops)
//#define NOROWPTRS
//#define ATYPE double *
//#define AROW(i) (m_A+(i)*m_nskip)

//***************************************************************************

#define dMIN(A,B)  ((A)>(B) ? (B) : (A))
#define dMAX(A,B)  ((B)>(A) ? (B) : (A))

#define LMATRIX_ALIGNMENT       dMAX(64, EFFICIENT_ALIGNMENT)

//***************************************************************************

// transfer b-values to x-values
template<bool zero_b>
inline 
void transfer_b_to_x(double pairsbx[PBX__MAX], unsigned n)
{
    double *const endbx = pairsbx + (size_t)n * PBX__MAX;
    for (double *currbx = pairsbx; currbx != endbx; currbx += PBX__MAX) {
        currbx[PBX_X] = currbx[PBX_B];
        if (zero_b) {
            currbx[PBX_B] = 0.0;
        }
    }
}

// swap row/column i1 with i2 in the n*n matrix A. the leading dimension of
// A is nskip. this only references and swaps the lower triangle.
// if `do_fast_row_swaps' is nonzero and row pointers are being used, then
// rows will be swapped by exchanging row pointers. otherwise the data will
// be copied.

static 
void swapRowsAndCols (ATYPE A, unsigned n, unsigned i1, unsigned i2, unsigned nskip, 
                             int do_fast_row_swaps)
{
    assert (A && n > 0 && i1 >= 0 && i2 >= 0 && i1 < n && i2 < n &&
        nskip >= n && i1 < i2);

# ifdef ROWPTRS
    double *A_i1 = A[i1];
    double *A_i2 = A[i2];
    for (unsigned i=i1+1; i<i2; ++i) {
        double *A_i_i1 = A[i] + i1;
        A_i1[i] = *A_i_i1;
        *A_i_i1 = A_i2[i];
    }
    A_i1[i2] = A_i1[i1];
    A_i1[i1] = A_i2[i1];
    A_i2[i1] = A_i2[i2];
    // swap rows, by swapping row pointers
    if (do_fast_row_swaps) {
        A[i1] = A_i2;
        A[i2] = A_i1;
    }
    else {
        // Only swap till i2 column to match A plain storage variant.
        for (unsigned k = 0; k <= i2; ++k) {
            std::swap(A_i1[k], A_i2[k]);
        }
    }
    // swap columns the hard way
    for (unsigned j = i2 + 1; j < n; ++j) {
        double *A_j = A[j];
        std::swap(A_j[i1], A_j[i2]);
    }
# else
    double *A_i1 = A + (size_t)nskip * i1;
    double *A_i2 = A + (size_t)nskip * i2;

    for (unsigned k = 0; k < i1; ++k) {
        std::swap(A_i1[k], A_i2[k]);
    }

    double *A_i = A_i1 + nskip;
    for (unsigned i= i1 + 1; i < i2; A_i += nskip, ++i) {
        std::swap(A_i2[i], A_i[i1]);
    }

    std::swap(A_i1[i1], A_i2[i2]);

    double *A_j = A_i2 + nskip;
    for (unsigned j = i2 + 1; j < n; A_j += nskip, ++j) {
        std::swap(A_j[i1], A_j[i2]);
    }
# endif
}


// swap two indexes in the n*n LCP problem. i1 must be <= i2.

static 
void swapProblem (ATYPE A, double pairsbx[PBX__MAX], double *w, double pairslh[PLH__MAX],
                         unsigned *p, bool *state, int *findex,
                         unsigned n, unsigned i1, unsigned i2, unsigned nskip,
                         int do_fast_row_swaps)
{
    assert (n>0 && i1 < n && i2 < n && nskip >= n && i1 <= i2);
    
    if (i1 != i2) {
        swapRowsAndCols (A, n, i1, i2, nskip, do_fast_row_swaps);

        std::swap((pairsbx + (size_t)i1 * PBX__MAX)[PBX_B], (pairsbx + (size_t)i2 * PBX__MAX)[PBX_B]);
        std::swap((pairsbx + (size_t)i1 * PBX__MAX)[PBX_X], (pairsbx + (size_t)i2 * PBX__MAX)[PBX_X]);
        assert(PBX__MAX == 2);

        std::swap(w[i1], w[i2]);

        std::swap((pairslh + (size_t)i1 * PLH__MAX)[PLH_LO], (pairslh + (size_t)i2 * PLH__MAX)[PLH_LO]);
        std::swap((pairslh + (size_t)i1 * PLH__MAX)[PLH_HI], (pairslh + (size_t)i2 * PLH__MAX)[PLH_HI]);
        assert(PLH__MAX == 2);

        std::swap(p[i1], p[i2]);
        std::swap(state[i1], state[i2]);

        if (findex != NULL) {
            std::swap(findex[i1], findex[i2]);
        }
    }
}

//***************************************************************************
// dLCP manipulator object. this represents an n*n LCP problem.
//
// two index sets C and N are kept. each set holds a subset of
// the variable indexes 0..n-1. an index can only be in one set.
// initially both sets are empty.
//
// the index set C is special: solutions to A(C,C)\A(C,i) can be generated.

//***************************************************************************
// fast implementation of dLCP. see the above definition of dLCP for
// interface comments.
//
// `p' records the permutation of A,x,b,w,etc. p is initially 1:n and is
// permuted as the other vectors/matrices are permuted.
//
// A,x,b,w,lo,hi,state,findex,p,c are permuted such that sets C,N have
// contiguous indexes. the don't-care indexes follow N.
//
// an L*D*L' factorization is maintained of A(C,C), and whenever indexes are
// added or removed from the set C the factorization is updated.
// thus L*D*L'=A[C,C], i.e. a permuted top left nC*nC submatrix of A.
// the leading dimension of the matrix L is always `nskip'.
//
// at the start there may be other indexes that are unbounded but are not
// included in `nub'. dLCP will permute the matrix so that absolutely all
// unbounded vectors are at the start. thus there may be some initial
// permutation.
//
// the algorithms here assume certain patterns, particularly with respect to
// index transfer.

#ifdef dLCP_FAST

struct dLCP {
    const unsigned m_n;
    const unsigned m_nskip;
    unsigned m_nub;
    unsigned m_nC, m_nN;				// size of each index set
    ATYPE const m_A;				// A rows
    double *const m_pairsbx, *const m_w, *const m_pairslh;	// permuted LCP problem data
    double *const m_L, *const m_d;				// L*D*L' factorization of set C
    double *const m_Dell, *const m_ell, *const m_tmp;
    bool *const m_state;
    int *const m_findex;
    unsigned *const m_p, *const m_C;

    dLCP (unsigned _n, unsigned _nskip, unsigned _nub, double *_Adata, double *_pairsbx, double *_w,
        double *_pairslh, double *_L, double *_d,
        double *_Dell, double *_ell, double *_tmp,
        bool *_state, int *_findex, unsigned *_p, unsigned *_C, double **Arows);
    unsigned getNub() const { return m_nub; }
    void transfer_i_to_C (unsigned i);
    void transfer_i_to_N (unsigned /*i*/) { m_nN++; }			// because we can assume C and N span 1:i-1
    void transfer_i_from_N_to_C (unsigned i);
    void transfer_i_from_C_to_N (unsigned i, void *tmpbuf);
    static size_t estimate_transfer_i_from_C_to_N_mem_req(unsigned nC, unsigned nskip) { return dEstimateLDLTRemoveTmpbufSize(nC, nskip); }
    unsigned numC() const { return m_nC; }
    unsigned numN() const { return m_nN; }
    unsigned indexC (unsigned i) const { return i; }
    unsigned indexN (unsigned i) const { return i+m_nC; }
    double Aii (unsigned i) const  { return AROW(i)[i]; }
    template<unsigned q_stride>
    double AiC_times_qC (unsigned i, double *q) const { return calculateLargeVectorDot<q_stride> (AROW(i), q, m_nC); }
    template<unsigned q_stride>
    double AiN_times_qN (unsigned i, double *q) const { return calculateLargeVectorDot<q_stride> (AROW(i) + m_nC, q + (size_t)m_nC * q_stride, m_nN); }
    void pN_equals_ANC_times_qC (double *p, double *q);
    void pN_plusequals_ANi (double *p, unsigned i, bool dir_positive);
    template<unsigned p_stride>
    void pC_plusequals_s_times_qC (double *p, double s, double *q);
    void pN_plusequals_s_times_qN (double *p, double s, double *q);
    void solve1 (double *a, unsigned i, bool dir_positive, int only_transfer=0);
    void unpermute_X();
    void unpermute_W();
};


dLCP::dLCP (unsigned _n, unsigned _nskip, unsigned _nub, double *_Adata, double *_pairsbx, double *_w,
            double *_pairslh, double *_L, double *_d,
            double *_Dell, double *_ell, double *_tmp,
            bool *_state, int *_findex, unsigned *_p, unsigned *_C, double **Arows):
    m_n(_n), m_nskip(_nskip), m_nub(_nub), m_nC(0), m_nN(0),
# ifdef ROWPTRS
    m_A(Arows),
#else
    m_A(_Adata),
#endif
    m_pairsbx(_pairsbx), m_w(_w), m_pairslh(_pairslh), 
    m_L(_L), m_d(_d), m_Dell(_Dell), m_ell(_ell), m_tmp(_tmp),
    m_state(_state), m_findex(_findex), m_p(_p), m_C(_C)
{
    dxtSetZero<PBX__MAX>(m_pairsbx + PBX_X, m_n);

    {
# ifdef ROWPTRS
        // make matrix row pointers
        double *aptr = _Adata;
        ATYPE A = m_A;
        const unsigned n = m_n, nskip = m_nskip;
        for (unsigned k=0; k<n; aptr+=nskip, ++k) A[k] = aptr;
# endif
    }

    {
        unsigned *p = m_p;
        const unsigned n = m_n;
        for (unsigned k=0; k != n; ++k) p[k] = k;		// initially unpermutted
    }

    /*
    // for testing, we can do some random swaps in the area i > nub
    {
    const unsigned n = m_n;
    const unsigned nub = m_nub;
    if (nub < n) {
    for (unsigned k=0; k<100; k++) {
    unsigned i1,i2;
    do {
    i1 = dRandInt(n-nub)+nub;
    i2 = dRandInt(n-nub)+nub;
    }
    while (i1 > i2); 
    //printf ("--> %d %d\n",i1,i2);
    swapProblem (m_A, m_pairsbx, m_w, m_pairslh, m_p, m_state, m_findex, n, i1, i2, m_nskip, 0);
    }
    }
    */

    // permute the problem so that *all* the unbounded variables are at the
    // start, i.e. look for unbounded variables not included in `nub'. we can
    // potentially push up `nub' this way and get a bigger initial factorization.
    // note that when we swap rows/cols here we must not just swap row pointers,
    // as the initial factorization relies on the data being all in one chunk.
    // variables that have findex >= 0 are *not* considered to be unbounded even
    // if lo=-inf and hi=inf - this is because these limits may change during the
    // solution process.

    {
        int *findex = m_findex;
        double *pairslh = m_pairslh;
        const unsigned n = m_n;
        for (unsigned k = m_nub; k < n; ++k) {
            if (findex && findex[k] >= 0) continue;
            if ((pairslh + (size_t)k * PLH__MAX)[PLH_LO] == -INFINITY && (pairslh + (size_t)k * PLH__MAX)[PLH_HI] == INFINITY) {
                swapProblem (m_A, m_pairsbx, m_w, pairslh, m_p, m_state, findex, n, m_nub, k, m_nskip, 0);
                m_nub++;
            }
        }
    }

    // if there are unbounded variables at the start, factorize A up to that
    // point and solve for x. this puts all indexes 0..nub-1 into C.
    if (m_nub > 0) {
        const unsigned nub = m_nub;
        {
            double *Lrow = m_L;
            const unsigned nskip = m_nskip;
            for (unsigned j = 0; j < nub; Lrow += nskip, ++j) memcpy(Lrow, AROW(j), (j + 1) * sizeof(double));
        }
        transfer_b_to_x<false> (m_pairsbx, nub);
        factorMatrixAsLDLT<1> (m_L, m_d, nub, m_nskip);
        solveEquationSystemWithLDLT<1, PBX__MAX> (m_L, m_d, m_pairsbx + PBX_X, nub, m_nskip);
        dSetZero (m_w, nub);
        {
            unsigned *C = m_C;
            for (unsigned k = 0; k < nub; ++k) C[k] = k;
        }
        m_nC = nub;
    }

    // permute the indexes > nub such that all findex variables are at the end
    if (m_findex) {
        const unsigned nub = m_nub;
        int *findex = m_findex;
        unsigned num_at_end = 0;
        for (unsigned k = m_n; k > nub; ) {
            --k;
            if (findex[k] >= 0) {
                swapProblem (m_A, m_pairsbx, m_w, m_pairslh, m_p, m_state, findex, m_n, k, m_n - 1 - num_at_end, m_nskip, 1);
                num_at_end++;
            }
        }
    }

    // print info about indexes
    /*
    {
    const unsigned n = m_n;
    const unsigned nub = m_nub;
    for (unsigned k=0; k<n; k++) {
    if (k<nub) printf ("C");
    else if ((m_pairslh + (size_t)k * PLH__MAX)[PLH_LO] == -dInfinity && (m_pairslh + (size_t)k * PLH__MAX)[PLH_HI] == dInfinity) printf ("c");
    else printf (".");
    }
    printf ("\n");
    }
    */
}

void dLCP::transfer_i_to_C (unsigned i)
{
    {
        const unsigned nC = m_nC;

        if (nC > 0) {
            // ell,Dell were computed by solve1(). note, ell = D \ L1solve (L,A(i,C))
            double *const Ltgt = m_L + (size_t)m_nskip * nC, *ell = m_ell;
            memcpy(Ltgt, ell, nC * sizeof(double));

            double ell_Dell_dot = dxDot(m_ell, m_Dell, nC);
            double AROW_i_i = AROW(i)[i] != ell_Dell_dot ? AROW(i)[i] : std::nextafter(AROW(i)[i], INFINITY); // A hack to avoid getting a zero in the denominator
            m_d[nC] = 1.0 / (AROW_i_i - ell_Dell_dot);
        }
        else {
            m_d[0] = 1.0 / (AROW(i)[i]);
        }

        swapProblem (m_A, m_pairsbx, m_w, m_pairslh, m_p, m_state, m_findex, m_n, nC, i, m_nskip, 1);

        m_C[nC] = nC;
        m_nC = nC + 1; // nC value is outdated after this line
    }

}


void dLCP::transfer_i_from_N_to_C (unsigned i)
{
    {
        const unsigned nC = m_nC;
        if (nC > 0) {
            {
                double *const aptr = AROW(i);
                double *Dell = m_Dell;
                const unsigned *C = m_C;
#   ifdef NUB_OPTIMIZATIONS
                // if nub>0, initial part of aptr unpermuted
                const unsigned nub = m_nub;
                unsigned j=0;
                for ( ; j<nub; ++j) Dell[j] = aptr[j];
                for ( ; j<nC; ++j) Dell[j] = aptr[C[j]];
#   else
                for (unsigned j=0; j<nC; ++j) Dell[j] = aptr[C[j]];
#   endif
            }
            solveL1Straight<1>(m_L, m_Dell, nC, m_nskip);

            double ell_Dell_dot = 0.0;
            double *const Ltgt = m_L + (size_t)m_nskip * nC;
            double *ell = m_ell, *Dell = m_Dell, *d = m_d;
            for (unsigned j = 0; j < nC; ++j) {
                double ell_j, Dell_j = Dell[j];
                Ltgt[j] = ell[j] = ell_j = Dell_j * d[j];
                ell_Dell_dot += ell_j * Dell_j;
            }
            
            double AROW_i_i = AROW(i)[i] != ell_Dell_dot ? AROW(i)[i] : std::nextafter(AROW(i)[i], INFINITY); // A hack to avoid getting a zero in the denominator
            m_d[nC] = 1.0 / (AROW_i_i - ell_Dell_dot);
        }
        else {
            m_d[0] = 1.0 / (AROW(i)[i]);
        }

        swapProblem (m_A, m_pairsbx, m_w, m_pairslh, m_p, m_state, m_findex, m_n, nC, i, m_nskip, 1);

        m_C[nC] = nC;
        m_nN--;
        m_nC = nC + 1; // nC value is outdated after this line
    }

    // @@@ TO DO LATER
    // if we just finish here then we'll go back and re-solve for
    // delta_x. but actually we can be more efficient and incrementally
    // update delta_x here. but if we do this, we wont have ell and Dell
    // to use in updating the factorization later.

}


void dLCP::transfer_i_from_C_to_N (unsigned i, void *tmpbuf)
{
    {
        unsigned *C = m_C;
        // remove a row/column from the factorization, and adjust the
        // indexes (black magic!)
        int last_idx = -1;
        const unsigned nC = m_nC;
        unsigned j = 0;
        for ( ; j < nC; ++j) {
            if (C[j] == nC - 1) {
                last_idx = j;
            }
            if (C[j] == i) {
                dxLDLTRemove (m_A, C, m_L, m_d, m_n, nC, j, m_nskip, tmpbuf);
                unsigned k;
                if (last_idx == -1) {
                    for (k = j + 1 ; k < nC; ++k) {
                        if (C[k] == nC - 1) {
                            break;
                        }
                    }
                    assert (k < nC);
                }
                else {
                    k = last_idx;
                }
                C[k] = C[j];
                if (j != (nC - 1)) memmove (C + j, C + j + 1, (nC - j - 1) * sizeof(C[0]));
                break;
            }
        }
        assert (j < nC);

        swapProblem (m_A, m_pairsbx, m_w, m_pairslh, m_p, m_state, m_findex, m_n, i, nC - 1, m_nskip, 1);

        m_nN++;
        m_nC = nC - 1; // nC value is outdated after this line
    }

}


void dLCP::pN_equals_ANC_times_qC (double *p, double *q)
{
    // we could try to make this matrix-vector multiplication faster using
    // outer product matrix tricks, e.g. with the dMultidotX() functions.
    // but i tried it and it actually made things slower on random 100x100
    // problems because of the overhead involved. so we'll stick with the
    // simple method for now.
    const unsigned nC = m_nC;
    double *ptgt = p + nC;
    const unsigned nN = m_nN;
    for (unsigned i = 0; i < nN; ++i) {
        ptgt[i] = dxDot (AROW(i + nC), q, nC);
    }
}

void dLCP::pN_plusequals_ANi (double *p, unsigned i, bool dir_positive)
{
    const unsigned nC = m_nC;
    double *aptr = AROW(i) + nC;
    double *ptgt = p + nC;
    if (dir_positive) {
        const unsigned nN = m_nN;
        for (unsigned j=0; j < nN; ++j) ptgt[j] += aptr[j];
    }
    else {
        const unsigned nN = m_nN;
        for (unsigned j=0; j < nN; ++j) ptgt[j] -= aptr[j];
    }
}

template<unsigned p_stride>
void dLCP::pC_plusequals_s_times_qC (double *p, double s, double *q)
{
    const unsigned nC = m_nC;
    double *q_end = q + nC;
    for (; q != q_end; p += p_stride, ++q) {
        *p += s * (*q);
    }
}

void dLCP::pN_plusequals_s_times_qN (double *p, double s, double *q)
{
    const unsigned nC = m_nC;
    double *ptgt = p + nC, *qsrc = q + nC;
    const unsigned nN = m_nN;
    for (unsigned i = 0; i < nN; ++i) {
        ptgt[i] += s * qsrc[i];
    }
}

void dLCP::solve1 (double *a, unsigned i, bool dir_positive, int only_transfer)
{
    // the `Dell' and `ell' that are computed here are saved. if index i is
    // later added to the factorization then they can be reused.
    //
    // @@@ question: do we need to solve for entire delta_x??? yes, but
    //     only if an x goes below 0 during the step.

    const unsigned nC = m_nC;
    if (nC > 0) {
        {
            double *Dell = m_Dell;
            unsigned *C = m_C;
            double *aptr = AROW(i);
#   ifdef NUB_OPTIMIZATIONS
            // if nub>0, initial part of aptr[] is guaranteed unpermuted
            const unsigned nub = m_nub;
            unsigned j = 0;
            for ( ; j < nub; ++j) Dell[j] = aptr[j];
            for ( ; j < nC; ++j) Dell[j] = aptr[C[j]];
#   else
            for (unsigned j = 0; j < nC; ++j) Dell[j] = aptr[C[j]];
#   endif
        }
        solveL1Straight<1>(m_L, m_Dell, nC, m_nskip);
        {
            double *ell = m_ell, *Dell = m_Dell, *d = m_d;
            for (unsigned j = 0; j < nC; ++j) ell[j] = Dell[j] * d[j];
        }

        if (!only_transfer) {
            double *tmp = m_tmp, *ell = m_ell;
            {
                for (unsigned j = 0; j < nC; ++j) tmp[j] = ell[j];
            }
            solveL1Transposed<1>(m_L, tmp, nC, m_nskip);
            if (dir_positive) {
                unsigned *C = m_C;
                double *tmp = m_tmp;
                for (unsigned j = 0; j < nC; ++j) a[C[j]] = -tmp[j];
            } else {
                unsigned *C = m_C;
                double *tmp = m_tmp;
                for (unsigned j = 0; j < nC; ++j) a[C[j]] = tmp[j];
            }
        }
    }
}


void dLCP::unpermute_X()
{
    unsigned *p = m_p;
    double *pairsbx = m_pairsbx;
    const unsigned n = m_n;
    for (unsigned j = 0; j < n; ++j) {
        unsigned k = p[j];
        if (k != j) {
            // p[j] = j; -- not going to be checked anymore anyway
            double x_j = (pairsbx + (size_t)j * PBX__MAX)[PBX_X];
            for (;;) {
                std::swap(x_j, (pairsbx + (size_t)k * PBX__MAX)[PBX_X]);

                unsigned orig_k = p[k];
                p[k] = k;
                if (orig_k == j) {
                    break;
                }
                k = orig_k;
            }
            (pairsbx + (size_t)j * PBX__MAX)[PBX_X] = x_j;
        }
    }
}

void dLCP::unpermute_W()
{
    memcpy (m_tmp, m_w, m_n * sizeof(double));

    const unsigned *p = m_p;
    double *w = m_w, *tmp = m_tmp;
    const unsigned n = m_n;
    for (unsigned j = 0; j < n; ++j) {
        unsigned k = p[j];
        w[k] = tmp[j];
    }
}

#endif // dLCP_FAST


//***************************************************************************
// if all the variables are unbounded then we can just factor, solve, and return

void dxSolveLCP_AllUnbounded (LCPMemoryManager* lcpMemoryManager, int n) {
    int nskip = lcpMemoryManager->nPad(n);
    double *A = lcpMemoryManager->A.get_ref();
    double *bx = lcpMemoryManager->bx.get_ref();

/*
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%lf\t", A[i * nskip + j]);
        }
        printf("\n");
    }

    for (int i = 0; i < n; i++)
        printf("%lf\n", pairsbx[2 * i + PBX_B]);
*/

    transfer_b_to_x<true>(bx, n);

    factorMatrixAsLDLT<PBX__MAX> (A, bx + PBX_B, n, nskip);
    solveEquationSystemWithLDLT<PBX__MAX, PBX__MAX> (A, bx + PBX_B, bx + PBX_X, n, nskip);
}

//***************************************************************************
// an optimized Dantzig LCP driver routine for the lo-hi LCP problem.

static 
void dxSolveLCP_Generic (LCPMemoryManager* lcpMemoryManager, int n, int nub) {
    int nskip = lcpMemoryManager->nPad(n);
    double *A = lcpMemoryManager->A.get_ref();
    double *pairsbx = lcpMemoryManager->bx.get_ref();
    double *pairslh = lcpMemoryManager->lohi.get_ref();
    int *findex = lcpMemoryManager->findex.get_ref();

    assert (n > 0 && A && pairsbx && pairslh && nub >= 0 && nub < n);
# ifndef dNODEBUG
    {
        // check restrictions on lo and hi
        double *endlh = pairslh + (size_t)n * PLH__MAX;
//        for (double *currlh = pairslh; currlh != endlh; currlh += PLH__MAX) printf ("%lf %lf\n", currlh[PLH_LO], currlh[PLH_HI]);

        for (double *currlh = pairslh; currlh != endlh; currlh += PLH__MAX) assert (currlh[PLH_LO] <= 0 && currlh[PLH_HI] >= 0);
    }
# endif

    double *L = lcpMemoryManager->L.get_ref();
    double *d = lcpMemoryManager->d.get_ref();
    double *w = lcpMemoryManager->w.get_ref();
    double *delta_w = lcpMemoryManager->delta_w.get_ref();
    double *delta_x = lcpMemoryManager->delta_x.get_ref();
    double *Dell = lcpMemoryManager->Dell.get_ref();
    double *ell = lcpMemoryManager->ell.get_ref();
#ifdef ROWPTRS
    double **Arows = lcpMemoryManager->Arows.get_ref();
#else
    double **Arows = NULL;
#endif
    unsigned *p = lcpMemoryManager->p.get_ref();
    unsigned *C = lcpMemoryManager->C.get_ref();

    // for i in N, state[i] is 0 if x(i)==lo(i) or 1 if x(i)==hi(i)
    bool *state = lcpMemoryManager->state.get_ref();

    // create LCP object. note that tmp is set to delta_w to save space, this
    // optimization relies on knowledge of how tmp is used, so be careful!
    dLCP lcp(n, nskip, nub, A, pairsbx, w, pairslh, L, d, Dell, ell, delta_w, state, findex, p, C, Arows);
    unsigned adj_nub = lcp.getNub();

    // loop over all indexes adj_nub..n-1. for index i, if x(i),w(i) satisfy the
    // LCP conditions then i is added to the appropriate index set. otherwise
    // x(i),w(i) is driven either +ve or -ve to force it to the valid region.
    // as we drive x(i), x(C) is also adjusted to keep w(C) at zero.
    // while driving x(i) we maintain the LCP conditions on the other variables
    // 0..i-1. we do this by watching out for other x(i),w(i) values going
    // outside the valid region, and then switching them between index sets
    // when that happens.

    bool hit_first_friction_index = false;
    for (unsigned i = adj_nub; i < n; ++i) {
        bool s_error = false;
        // the index i is the driving index and indexes i+1..n-1 are "dont care",
        // i.e. when we make changes to the system those x's will be zero and we
        // don't care what happens to those w's. in other words, we only consider
        // an (i+1)*(i+1) sub-problem of A*x=b+w.

        // if we've hit the first friction index, we have to compute the lo and
        // hi values based on the values of x already computed. we have been
        // permuting the indexes, so the values stored in the findex vector are
        // no longer valid. thus we have to temporarily unpermute the x vector. 
        // for the purposes of this computation, 0*infinity = 0 ... so if the
        // contact constraint's normal force is 0, there should be no tangential
        // force applied.

        if (!hit_first_friction_index && findex && findex[i] >= 0) {
            // un-permute x into delta_w, which is not being used at the moment
            for (unsigned j = 0; j < n; ++j) delta_w[p[j]] = (pairsbx + (size_t)j * PBX__MAX)[PBX_X];

            // set lo and hi values
            for (unsigned k = i; k < n; ++k) {
                double *currlh = pairslh + (size_t)k * PLH__MAX;
                double wfk = delta_w[findex[k]];
                if (wfk == 0) {
                    currlh[PLH_HI] = 0;
                    currlh[PLH_LO] = 0;
                }
                else {
                    currlh[PLH_HI] = fabs (currlh[PLH_HI] * wfk);
                    currlh[PLH_LO] = -currlh[PLH_HI];
                }
            }
            hit_first_friction_index = true;
        }

        // thus far we have not even been computing the w values for indexes
        // greater than i, so compute w[i] now.
        double wPrep = lcp.AiC_times_qC<PBX__MAX> (i, pairsbx + PBX_X) + lcp.AiN_times_qN<PBX__MAX> (i, pairsbx + PBX_X);

        double *currbx = pairsbx + (size_t)i * PBX__MAX;

        w[i] = wPrep - currbx[PBX_B];

        // if lo=hi=0 (which can happen for tangential friction when normals are
        // 0) then the index will be assigned to set N with some state. however,
        // set C's line has zero size, so the index will always remain in set N.
        // with the "normal" switching logic, if w changed sign then the index
        // would have to switch to set C and then back to set N with an inverted
        // state. this is pointless, and also computationally expensive. to
        // prevent this from happening, we use the rule that indexes with lo=hi=0
        // will never be checked for set changes. this means that the state for
        // these indexes may be incorrect, but that doesn't matter.

        double *currlh = pairslh + (size_t)i * PLH__MAX;

        // see if x(i),w(i) is in a valid region
        if (currlh[PLH_LO] == 0 && w[i] >= 0) {
            lcp.transfer_i_to_N (i);
            state[i] = false;
        }
        else if (currlh[PLH_HI] == 0 && w[i] <= 0) {
            lcp.transfer_i_to_N (i);
            state[i] = true;
        }
        else if (w[i] == 0) {
            // this is a degenerate case. by the time we get to this test we know
            // that lo != 0, which means that lo < 0 as lo is not allowed to be +ve,
            // and similarly that hi > 0. this means that the line segment
            // corresponding to set C is at least finite in extent, and we are on it.
            // NOTE: we must call lcp.solve1() before lcp.transfer_i_to_C()
            lcp.solve1 (delta_x, i, false, 1);

            lcp.transfer_i_to_C (i);
        }
        else {
            // we must push x(i) and w(i)
            for (;;) {
                // find direction to push on x(i)
                bool dir_positive = (w[i] <= 0);

                // compute: delta_x(C) = -dir*A(C,C)\A(C,i)
                lcp.solve1 (delta_x, i, dir_positive);

                // note that delta_x[i] = (dir_positive ? 1 : -1), but we wont bother to set it

                // compute: delta_w = A*delta_x ... note we only care about
                // delta_w(N) and delta_w(i), the rest is ignored
                lcp.pN_equals_ANC_times_qC (delta_w, delta_x);
                lcp.pN_plusequals_ANi (delta_w, i, dir_positive);
                delta_w[i] = dir_positive 
                    ? lcp.AiC_times_qC<1> (i, delta_x) + lcp.Aii(i)
                    : lcp.AiC_times_qC<1> (i, delta_x) - lcp.Aii(i);

                // find largest step we can take (size=s), either to drive x(i),w(i)
                // to the valid LCP region or to drive an already-valid variable
                // outside the valid region.

                int cmd = 1;		// index switching command
                unsigned si = 0;		// si = index to switch if cmd>3

                double s = delta_w[i] != 0.0
                    ? -w[i] / delta_w[i]
                    : (w[i] != 0.0 ? copysign(INFINITY, -w[i]) : 0.0);
                    
                if (dir_positive) {
                    if (currlh[PLH_HI] < INFINITY) {
                        double s2 = (currlh[PLH_HI] - currbx[PBX_X]);	// was (hi[i]-x[i])/dirf	// step to x(i)=hi(i)
                        if (s2 < s) {
                            s = s2;
                            cmd = 3;
                        }
                    }
                }
                else {
                    if (currlh[PLH_LO] > -INFINITY) {
                        double s2 = (currbx[PBX_X] - currlh[PLH_LO]); // was (lo[i]-x[i])/dirf	// step to x(i)=lo(i)
                        if (s2 < s) {
                            s = s2;
                            cmd = 2;
                        }
                    }
                }

                {
                    const unsigned numN = lcp.numN();
                    for (unsigned k = 0; k < numN; ++k) {
                        const unsigned indexN_k = lcp.indexN(k);
                        if (!state[indexN_k] ? delta_w[indexN_k] < 0 : delta_w[indexN_k] > 0) {
                            // don't bother checking if lo=hi=0
                            double *indexlh = pairslh + (size_t)indexN_k * PLH__MAX;
                            if (indexlh[PLH_LO] == 0 && indexlh[PLH_HI] == 0) continue;
                            double s2 = -w[indexN_k] / delta_w[indexN_k];
                            if (s2 < s) {
                                s = s2;
                                cmd = 4;
                                si = indexN_k;
                            }
                        }
                    }
                }

                {
                    const unsigned numC = lcp.numC();
                    for (unsigned k = adj_nub; k < numC; ++k) {
                        const unsigned indexC_k = lcp.indexC(k);
                        double *indexlh = pairslh + (size_t)indexC_k * PLH__MAX;
                        if (delta_x[indexC_k] < 0 && indexlh[PLH_LO] > -INFINITY) {
                            double s2 = (indexlh[PLH_LO] - (pairsbx + (size_t)indexC_k * PBX__MAX)[PBX_X]) / delta_x[indexC_k];
                            if (s2 < s) {
                                s = s2;
                                cmd = 5;
                                si = indexC_k;
                            }
                        }
                        if (delta_x[indexC_k] > 0 && indexlh[PLH_HI] < INFINITY) {
                            double s2 = (indexlh[PLH_HI] - (pairsbx + (size_t)indexC_k * PBX__MAX)[PBX_X]) / delta_x[indexC_k];
                            if (s2 < s) {
                                s = s2;
                                cmd = 6;
                                si = indexC_k;
                            }
                        }
                    }
                }

                //static char* cmdstring[8] = {0,"->C","->NL","->NH","N->C",
                //			     "C->NL","C->NH"};
                //printf ("cmd=%d (%s), si=%d\n",cmd,cmdstring[cmd],(cmd>3) ? si : i);

                // if s <= 0 then we've got a problem. if we just keep going then
                // we're going to get stuck in an infinite loop. instead, just cross
                // our fingers and exit with the current solution.
                if (s <= 0.0) {
                    printf ("LCP internal error, s <= 0 (s=%.4e)",(double)s);
                    if (i < n) {
                        dxtSetZero<PBX__MAX>(currbx + PBX_X, n - i);
                        dxSetZero (w + i, n - i);
                    }
                    s_error = true;
                    break;
                }

                // apply x = x + s * delta_x
                lcp.pC_plusequals_s_times_qC<PBX__MAX> (pairsbx + PBX_X, s, delta_x);
                currbx[PBX_X] = dir_positive 
                    ? currbx[PBX_X] + s
                    : currbx[PBX_X] - s;

                // apply w = w + s * delta_w
                lcp.pN_plusequals_s_times_qN (w, s, delta_w);
                w[i] += s * delta_w[i];

                void *tmpbuf;
                // switch indexes between sets if necessary
                switch (cmd) {
                case 1:		// done
                    w[i] = 0;
                    lcp.transfer_i_to_C (i);
                    break;
                case 2:		// done
                    currbx[PBX_X] = currlh[PLH_LO];
                    state[i] = false;
                    lcp.transfer_i_to_N (i);
                    break;
                case 3:		// done
                    currbx[PBX_X] = currlh[PLH_HI];
                    state[i] = true;
                    lcp.transfer_i_to_N (i);
                    break;
                case 4:		// keep going
                    w[si] = 0;
                    lcp.transfer_i_from_N_to_C (si);
                    break;
                case 5:		// keep going
                    (pairsbx + (size_t)si * PBX__MAX)[PBX_X] = (pairslh + (size_t)si * PLH__MAX)[PLH_LO];
                    state[si] = false;
                    tmpbuf = lcpMemoryManager->tmpbuf.get_ref();
                    lcp.transfer_i_from_C_to_N (si, tmpbuf);
                    break;
                case 6:		// keep going
                    (pairsbx + (size_t)si * PBX__MAX)[PBX_X] = (pairslh + (size_t)si * PLH__MAX)[PLH_HI];
                    state[si] = true;
                    tmpbuf = lcpMemoryManager->tmpbuf.get_ref();
                    lcp.transfer_i_from_C_to_N (si, tmpbuf);
                    break;
                }

                if (cmd <= 3) break;
            } // for (;;)
        } // else

        if (s_error) {
            break;
        }
    } // for (unsigned i = adj_nub; i < n; ++i)

    lcp.unpermute_X(); // This destroys p[] and must be done last
}

void solveLCP (LCPMemoryManager* lcpMemoryManager, int n, int nub) {
    if (nub >= n)
        dxSolveLCP_AllUnbounded (lcpMemoryManager, n);
    else
        dxSolveLCP_Generic (lcpMemoryManager, n, nub);
}

//print out the feedback torques, and see if they comply to the limits - and see what happens when torque limits are not active...