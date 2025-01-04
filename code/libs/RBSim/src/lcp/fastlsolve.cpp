#include "RBSim/lcp/fastldltsolve_impl.h"

void solveL1(const double *L, double *B, int n, int lskip1) {
    assert(n != 0);

    if (n != 0) {
        assert(L != NULL);
        assert(B != NULL);

        solveL1Straight<1>(L, B, n, lskip1);
    }
}

