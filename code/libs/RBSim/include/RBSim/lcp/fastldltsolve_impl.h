#pragma once

#include <assert.h>

#include "RBSim/lcp/fastlsolve_impl.h"
#include "RBSim/lcp/fastltsolve_impl.h"
#include "RBSim/lcp/fastvecscale_impl.h"


template<unsigned int d_stride, unsigned int b_stride>
void solveEquationSystemWithLDLT(const double *L, const double *d, double *b, unsigned rowCount, unsigned rowSkip)
{
    assert(L != NULL);
    assert(d != NULL);
    assert(b != NULL);
    assert(rowCount > 0);
    assert(rowSkip >= rowCount);

    solveL1Straight<b_stride>(L, b, rowCount, rowSkip);
    scaleLargeVector<b_stride, d_stride>(b, d, rowCount);
    solveL1Transposed<b_stride>(L, b, rowCount, rowSkip);
}
