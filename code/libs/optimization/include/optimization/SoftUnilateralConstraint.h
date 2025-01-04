#pragma once

#include <math.h>
#include <optimization/ObjectiveFunction.h>

/** class used to model unilateral constraints of the type x > l using a C2
 * penalty energy f(x).
 * - l is the lower limit that x needs to be greater than
 * - if x > l, then the energy of the constraint, its gradient and hessian are
 * all 0 (i.e. inactive)
 * - epsilon is the value away from the limit (how much larger should x be
 * compared to l) after which f(x) = 0
 * - stiffness controls the rate at which f(x) increases if x < l
 */
class SoftLowerLimitConstraint {
private:
    double a1, b1, c1, a2, b2, c2, d2, epsilon;

public:
    double limit = 0;

public:
    SoftLowerLimitConstraint(double l, double stiffness, double epsilon);

    virtual ~SoftLowerLimitConstraint();

    //based on the value of the stiffness and eps parameters, update the coeffs
    //of the piece-wise cubic polynomial that we're working with here
    void setPolynomialCoeffs(double stiffness, double eps);

    void setLimit(double l);
    double getEpsilon();
    double getStiffness();
    double getLimit();

    // comptue f(x)
    double evaluate(double x);

    // compute df/dx
    double computeDerivative(double x);

    // compute ddf/dxdx
    double computeSecondDerivative(double x);
};

/**
 * class used to model unilateral constraints of the type x < u using a C2
 * penalty energy f(x).
 * - u is the upper limit that x needs to be less than
 * - if x < u, then the energy of the constraint, its gradient and hessian are
 * all 0 (i.e. inactive)
 * - epsilon is the value away from the limit (how much smaller should x be
 * compared to u) after which f(x) = 0
 * - stiffness controls the rate at which f(x) increases if x > u
 */
class SoftUpperLimitConstraint {
private:
    double a1, b1, c1, a2, b2, c2, d2, epsilon;
    double limit = 0;

public:
    SoftUpperLimitConstraint(double l, double stiffness, double epsilon);

    virtual ~SoftUpperLimitConstraint();

    //based on the value of the stiffness and eps parameters, update the coeffs
    //of the piece-wise cubic polynomial that we're working with here
    void setPolynomialCoeffs(double stiffness, double eps);

    void setLimit(double l);
    double getEpsilon();
    double getStiffness();
    double getLimit();

    // compute f(x)
    double evaluate(double x);

    // compute df/dx
    double computeDerivative(double x);

    // compute ddf/dxdx
    double computeSecondDerivative(double x);
};

class SoftBoundConstraint {
protected:
    SoftLowerLimitConstraint l;
    SoftUpperLimitConstraint u;

public:
    SoftBoundConstraint(double lLimit, double uLimit,
                        double stiffness = 1.0, double relEps = 0.1);
    virtual ~SoftBoundConstraint();

    // comptue f(x)
    double evaluate(double x);

    // compute df/dx
    double computeDerivative(double x);

    // compute ddf/dxdx
    double computeSecondDerivative(double x);

    double getEPS();
    double getStiffness();
    void setStiffnessAndRelEpsilon(double stiffness, double eps);
};
