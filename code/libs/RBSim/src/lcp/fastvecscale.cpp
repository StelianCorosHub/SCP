#include "RBSim/lcp/fastvecscale_impl.h"

void dScaleVector(double *a, const double *d, int n) {
    scaleLargeVector<1, 1>(a, d, n);
}

