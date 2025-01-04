#pragma once

#include <math.h>

class SoftEqualityConstraint {
private:
    double _rhs, _stiffness;

public:
    SoftEqualityConstraint(double rhs, double stiffness = 1.)
        : _rhs(rhs), _stiffness(stiffness) {}

    virtual ~SoftEqualityConstraint() {}

    void setRHS(double rhs) { _rhs = rhs; }

    void setStiffness(double stiffness) { _stiffness = stiffness; }

    // compute f(x)
    double evaluate(double x) const {
        return 0.5 * _stiffness * pow(x - _rhs, 2.);
    }

    // compute df/dx
    double computeDerivative(double x) const { return _stiffness * (x - _rhs); }

    // compute ddf/dxdx
    double computeSecondDerivative() const { return _stiffness; }
};
