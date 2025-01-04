#pragma once

#include <utils/mathUtils.h>

#include "optimization/ObjectiveFunction.h"
#include "optimization/Minimizer.h"

/*
    Basic idea: compute a (descent) search direction, run a line search along it until a smaller function value is found.
*/
class LineSearchMinimizer : public Minimizer {
public:
    LineSearchMinimizer(double solveResidual = 1e-5,
                        double lineSearchStartingStepSize = 1.0,
                        int maxLineSearchSteps = 10) {
        this->solveResidual = solveResidual;
        this->lineSearchStartingStepSize = lineSearchStartingStepSize;
        this->maxLineSearchSteps = maxLineSearchSteps;
    }

    bool minimize(ObjectiveFunction *function, Vector &x, int maxIterations = 100) const override {
        // this will be the search direction
        Vector dx(x.size());

        if (printOutput)
            std::cout << "Initial objective function value: "
                      << function->evaluate(x) << std::endl;

        for (int i = 0; i < maxIterations; i++) {
            function->prepareForOptimizationStep(x);
            dx = getSearchDirection(function, x);

            if (dx.norm() < solveResidual) return true;

            doLineSearch(function, dx, x);

            if (printOutput)
                std::cout
                    << "Objective function value after optimization step: "
                    << function->evaluate(x) << std::endl;
        }
        return false;
    }

    /**
     * this function needs to provide a search (i.e. descent!) direction dx evaluated at the current solution x.
     */
    virtual Vector getSearchDirection(ObjectiveFunction *function, const Vector &x) const = 0;

    /**
     * given the objective `function` and search direction `dx`, update the
     * candidate `x` via a bisection line search. If no better function value is found
     * given the parameters of the line search, then x is not modified.
     */
    virtual void doLineSearch(ObjectiveFunction *function, const Vector &dx, Vector &x) const {
        // line search
        double alpha = lineSearchStartingStepSize;  // initial step size
        // these will be candidate solutions that we try out
        Vector xStart(x);
        double initialFunctionValue = function->evaluate(xStart);
        int lSteps = MAX(maxLineSearchSteps, 0);

        for (int j = 0; j <= lSteps; j++) {
            // try a step of size `alpha` in the search (descent!) direction
            x = xStart + dx * alpha;

            // if the new function value is greater than the initial one, we've
            // gone uphill. Reduce alpha and try again...
            double f = function->evaluate(x);

            if (printOutput)
                printf(
                    "\t >>> Line search: alpha = %lf, fBefore = %lf, fAfter = "
                    "%lf\n",
                    alpha, initialFunctionValue, f);

            if (!std::isfinite(f) || f > initialFunctionValue)
                alpha /= 2.0;
            else
                return;
        }
    }

public:
    double solveResidual = 1e-5;
    double lineSearchStartingStepSize = 1.0;
    int maxLineSearchSteps = 10;
};

class GradientDescentMinimizer : public LineSearchMinimizer {
public:
    GradientDescentMinimizer(double solveResidual = 1e-5,
                             double lineSearchStartingStepSize = 1.0,
                             int maxLineSearchSteps = 10)
        : LineSearchMinimizer(solveResidual, lineSearchStartingStepSize,
                              maxLineSearchSteps) {}

    virtual Vector getSearchDirection(ObjectiveFunction *function, const Vector &x) const {
        return function->getGradient(x) * -1;
    }
};
