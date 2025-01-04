#pragma once

#include "optimization/ObjectiveFunction.h"

class Minimizer {
public:
    /**
     * Returns true if a minimum of the objective `function` has been found. `x`
     * is the initial/current candidate, and will also store the next candidate
     * once the method has returned. The method will perform at most
     * maxIterations steps of the iterative optimization procedure.
     */
    virtual bool minimize(ObjectiveFunction *function, Vector &x, int maxIterations = 100) const = 0;

    bool printOutput = true;
};
