#pragma once

#include <utils/mathUtils.h>

#include <Eigen/Sparse>
#include <vector>

class ObjectiveFunction {
public:
    std::string description = "";

    ObjectiveFunction(std::string descr = "") { description = descr; }

    /**
     * this function is called before each and every optimization step starts --
     * can be used to precompute various quantities, set regularizers,
     * etc...
     */
    virtual void prepareForOptimizationStep(const Vector &x) {}

    /**
     * returns O(x)
     */
    virtual double evaluate(const Vector &x) = 0;

    virtual void computeGradient(const Vector &x, Vector &grad);

    Vector getGradient(const Vector &x) {
        Vector grad;
        computeGradient(x, grad);
        return grad;
    }

    virtual void computeHessian(const Vector &x, SparseMatrix &hessian);

    void addFiniteDifferenceGradientTo(const Vector &x, Vector &grad);

    void addFiniteDifferenceHessianEntriesTo(const Vector &x, std::vector<SMTriplet> &hessianEntries);

    void testGradientWithFD(const Vector &x);

    void testSparseHessianWithFD(const Vector &x, double maxAcceptableError = 10e-5, double maxAcceptableRelError = 10e-3);

//    void testSparseHessianPSD(const Vector &x);

//    void testDenseHessianWithFD(const Vector &x, double maxAcceptableError = 10e-5,
//                                 double maxAcceptableRelError = 10e-3);
//    void testDenseHessianPSD(const Vector &x);
//    virtual void computeHessian(const Vector &x, Matrix &hessian) const;
//    void addFiniteDifferenceHessianEntriesTo(const Vector &x, Matrix &hessian) const;

    virtual void printObjectiveDetails(const Vector &x) {}
};
