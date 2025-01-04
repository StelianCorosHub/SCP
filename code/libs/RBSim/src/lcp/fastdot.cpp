#include "RBSim/lcp/fastdot_impl.h"

double dxDot (const double *a, const double *b, unsigned n) {
    return calculateLargeVectorDot<1>(a, b, n);
}

double dDot (const double *a, const double *b, int n) {
    return dxDot (a, b, n);
}

