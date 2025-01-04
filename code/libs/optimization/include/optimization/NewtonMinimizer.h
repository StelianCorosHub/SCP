#pragma once

#include "optimization/GradientDescentMinimizer.h"
#include "optimization/ObjectiveFunction.h"

class NewtonMinimizer : public LineSearchMinimizer {
public:
    NewtonMinimizer(double solveResidual = 1e-5,
                    double lineSearchStartingStepSize = 1.0,
                    int maxLineSearchSteps = 10)
        : LineSearchMinimizer(solveResidual, lineSearchStartingStepSize,
                                 maxLineSearchSteps) {}

    virtual ~NewtonMinimizer() {}


    /**
     * With Newton's method, the search direction is given by - H^-1 * g
     */
    Vector getSearchDirection(ObjectiveFunction *function,
                                const Vector &x) const override {
        int maxSteps = MAX(maxDynamicRegularizationSteps, 0);
        double currentReg = reg;

        Vector dx(x.size());

        Vector r(x.size());

        // prepare gradient
        Vector g = function->getGradient(x);
        g *= -1;
        // and prepare hessian
        function->computeHessian(x, hessian);

        //		Matrix tmpH =
        // hessian.triangularView<Eigen::Lower>().toDense() +
        // hessian.transpose().triangularView<Eigen::StrictlyUpper>().toDense();
        //		Eigen::EigenSolver<Matrix> es(tmpH);
        //		std::cout << "The eigenvalues of A are:" << std::endl <<
        // es.eigenvalues() << std::endl; 		print("g.m", g);
        // print("H.m", hessian);

        for (int i = 0; i <= maxSteps; i++) {
            // add regularization
            if (currentReg > 0) {
                r.setConstant(currentReg);
                hessian += r.asDiagonal();
            }

            // solve for search direction: dx = -Hes^-1 * grad
            Eigen::SimplicialLDLT<SparseMatrix, Eigen::Lower> solver;
            solver.compute(hessian);
            dx = solver.solve(g);

//            if (printOutput) checkSymmetricLinearSystemResidual(hessian, dx, g);

            // check if we have a descent direction...
            double dotProduct = dx.dot(g);

            if (dotProduct <= 0) {
                // if not, increase regularizer and try again...
                if (printOutput) {
                    std::cout << "Search direction is NOT a descent direction";
                    std::cout << " (g.dx = " << dotProduct << ", ang = "
                              << DEG(safeACOS(dx.dot(g) /
                                              (dx.norm() * g.norm() + 1e-10)))
                              << " deg, reg = " << currentReg << ")"
                              << std::endl;
                }
                currentReg += 1e-4;
                currentReg *= 10;
            } else {
                // if it's a descent direction, we're good to go...
                if (printOutput) {
                    std::cout << "Search direction is a descent direction";
                    std::cout << " (g.dx = " << dotProduct << ", ang = "
                              << DEG(safeACOS(dx.dot(g) /
                                              (dx.norm() * g.norm() + 1e-10)))
                              << " deg, reg = " << currentReg << ")"
                              << std::endl;
                }
                break;
            }
        }

        return dx;
    }

public:
    mutable SparseMatrix hessian;
    double reg = 1e-4;
    int maxDynamicRegularizationSteps = 5;
};


