#include "RBSim/lcp/fastldltsolve_impl.h"

void solveLDLT(const double *L, const double *d, double *b, int n, int nskip) {
    assert(n != 0);

    if (n != 0) {
        assert(L != NULL);
        assert(d != NULL);
        assert(b != NULL);

        solveEquationSystemWithLDLT<1, 1>(L, d, b, n, nskip);
    }
}

