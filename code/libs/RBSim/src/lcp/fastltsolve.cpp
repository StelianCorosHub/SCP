#include "RBSim/lcp/fastldltsolve_impl.h"

void dSolveL1T(const double *L, double *B, int rowCount, int rowSkip) {
    assert(rowCount != 0);

    if (rowCount != 0) {
        assert(L != NULL);
        assert(B != NULL);

        solveL1Transposed<1>(L, B, rowCount, rowSkip);
    }
}
