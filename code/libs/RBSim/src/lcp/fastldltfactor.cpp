#include "RBSim/lcp/fastldltfactor_impl.h"

void factorLDLT(double *A, double *d, int n, int nskip1) {
    factorMatrixAsLDLT<1>(A, d, n, nskip1);
}
