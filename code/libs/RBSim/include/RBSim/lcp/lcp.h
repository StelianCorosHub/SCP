
#pragma once

#include <stddef.h>
#include <RBSim/lcp/memory.h>

enum dxLCPBXElement
{
    PBX__MIN,

    PBX_B = PBX__MIN,
    PBX_X,

    PBX__MAX,
};

enum dxLCPLHElement
{
    PLH__MIN,

    PLH_LO = PLH__MIN,
    PLH_HI,

    PLH__MAX,
};

template <typename T> class ManagedMemoryBlock {
private:
    T* mem = nullptr;
    int n = 0;

public:
    T* get_ref(int nNew) {
        alloc(nNew);
        return mem;
    }

    void alloc(int nNew) {
        if (n < nNew) {
            n = nNew;
            mem = (T*) realloc(mem, n * sizeof(T));
        }
    }

    T* get_ref() {
        return mem;
    }

    ManagedMemoryBlock() {}

    ~ManagedMemoryBlock() {
        delete mem;
    }

};

class LCPMemoryManager {
public:
    ManagedMemoryBlock<double> A;
    ManagedMemoryBlock<double> bx;
    ManagedMemoryBlock<double> lohi;
    ManagedMemoryBlock<int> findex;
    ManagedMemoryBlock<double> L;
    ManagedMemoryBlock<double> d;
    ManagedMemoryBlock<double> w;
    ManagedMemoryBlock<double> delta_w;
    ManagedMemoryBlock<double> delta_x;
    ManagedMemoryBlock<double> Dell;
    ManagedMemoryBlock<double> ell;
    ManagedMemoryBlock<double*> Arows;
    ManagedMemoryBlock<unsigned> p;
    ManagedMemoryBlock<unsigned> C;
    ManagedMemoryBlock<bool> state;
    ManagedMemoryBlock<double> tmpbuf;

public:
    LCPMemoryManager() {}

    ~LCPMemoryManager() {
    }

    // round an integer up to a multiple of 4, except that 0 and 1 are unmodified,
    // used to compute matrix leading dimensions
    double nPad(int a) {
        return dPAD(a);
    }

    void initialize(int n) {
        A.alloc(nPad(n) * n);
        bx.alloc(2 * n);
        lohi.alloc(2 * n);
        findex.alloc(n);
        d.alloc(n);
        w.alloc(n);
        delta_w.alloc(n);
        delta_x.alloc(n);
        Dell.alloc(n);
        ell.alloc(n);
        Arows.alloc(n);
        p.alloc(n);
        C.alloc(n);
        state.alloc(n);
        L.alloc(nPad(n) * n);
        tmpbuf.alloc(nPad(n) * n);
    }

};


void dxSolveLCP_AllUnbounded (LCPMemoryManager* lcpMemoryManager, int n);


/*
 * given (A,b,lo,hi), solve the LCP problem: A*x = b+w, where each x(i), w(i)
 * satisfies one of
    (1) x = lo, w >= 0
    (2) x = hi, w <= 0
    (3) lo < x < hi, w = 0

 * A is a matrix of dimension n*n, everything else is a vector of size n*1.
 * lo and hi can be +/- dInfinity as needed. the first `nub' variables are
 * unbounded, i.e. hi and lo are assumed to be +/- dInfinity.
 * we restrict lo(i) <= 0 and hi(i) >= 0.

 * If the `findex' (friction index) parameter is nonzero, it points to an array
 * of index values. in this case constraints that have findex[i] >= 0 are
 * special. all non-special constraints are solved for, then the lo and hi values
 * for the special constraints are set:
 *  hi[i] = abs( hi[i] * x[findex[i]] )
 *  lo[i] = -hi[i]
 *  and the solution continues. this mechanism allows a friction approximation
 *  to be implemented. the first `nub' variables are assumed to have findex < 0.
*/

void solveLCP (LCPMemoryManager* lcpMemoryManager, int n, int nub);
