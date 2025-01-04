#include "optimization/ObjectiveFunction.h"

#include "utils/logger.h"

void ObjectiveFunction::computeGradient(const Vector &x, Vector &grad) {
    resize(grad, (int)x.size());
    addFiniteDifferenceGradientTo(x, grad);
}

void ObjectiveFunction::computeHessian(const Vector &x,
                                   SparseMatrix &hessian) {
    hessian.resize((int)x.size(), (int)x.size());
    hessian.setZero();
    std::vector<SMTriplet> triplets;
    triplets.clear();
    addFiniteDifferenceHessianEntriesTo(x, triplets);
    hessian.setFromTriplets(triplets.begin(), triplets.end());
}

void ObjectiveFunction::addFiniteDifferenceGradientTo(const Vector &x,
                                                      Vector &grad) {
    const double h = 1e-5;
    for (int i = 0; i < x.size(); ++i) {
        Vector dx(x.size());
        dx.setZero();
        dx[i] = h;

        double Op = evaluate(x + dx);
        double Om = evaluate(x - dx);
        grad[i] += (Op - Om) / (2. * h);
    }
}

void ObjectiveFunction::addFiniteDifferenceHessianEntriesTo(
    const Vector &x, std::vector<SMTriplet> &hessianEntries) {
    const double h = 1e-7;
    for (int i = 0; i < x.size(); ++i) {
        Vector dx(x.size());
        dx.setZero();
        dx[i] = h;

        Vector gp(x.size()), gm(x.size());
        computeGradient(x + dx, gp);
        computeGradient(x - dx, gm);

        Vector hess = (gp - gm) / (2. * h);
        // only copy lower diagonal entries to the estimated hessian, as it is
        // supposed to be sparse...
        for (int j = i; j < x.size(); ++j) {
            if (abs(hess[j]) > 1e-12 || i == j)
                hessianEntries.push_back(SMTriplet(j, i, hess[j]));
        }
    }
}

void ObjectiveFunction::testGradientWithFD(const Vector &p) {
    double tol = 1e-4;
    double eps = 1e-10;

    Vector FDGradient;
    Vector analyticGradient;

    resize(FDGradient, (int)p.size());
    resize(analyticGradient, (int)p.size());

    addFiniteDifferenceGradientTo(p, FDGradient);
    computeGradient(p, analyticGradient);

    //	print("..\\out\\gradient_FD.m", FDGradient);
    //	print("..\\out\\gradient.m", analyticGradient);

    if (fabs(analyticGradient.norm() - FDGradient.norm()) > 0.00001) {
        Logger::print(
            "Objective Function: testing gradients...norms: "
            "analytic: %lf, FD: %lf\n",
            analyticGradient.norm(), FDGradient.norm());
    } else {
        Logger::print("Objective Function: testing gradients...\n");
    }

    for (int i = 0; i < p.size(); i++) {
        double absErr = std::abs(FDGradient[i] - analyticGradient[i]);
        double relError =
            2 * absErr /
            (eps + std::abs(analyticGradient[i]) + std::abs(FDGradient[i]));

        if (relError > tol && absErr > 1e-6) {
            Logger::print(Logger::RED,
                          "Mismatch element %d: Analytic val: %lf, FD val: "
                          "%lf. Error: %lf(%lf%%)\n",
                          i, analyticGradient[i], FDGradient[i], absErr,
                          relError * 100);
        }
    }
}

void ObjectiveFunction::testSparseHessianWithFD(const Vector &p,
                                          double maxAcceptableError,
                                          double maxAcceptableRelError) {
    double eps = 1e-10;

    SparseMatrix FDHessian(p.size(), p.size());
    SparseMatrix analyticHessian(p.size(), p.size());
    Array<SMTriplet> hessianEntries;

    addFiniteDifferenceHessianEntriesTo(p, hessianEntries);
    FDHessian.setFromTriplets(hessianEntries.begin(), hessianEntries.end());
    computeHessian(p, analyticHessian);
    Logger::print("Objective Function: testing hessians...\n");
    for (int i = 0; i < p.size(); i++) {
        for (
            int j = 0; j <= i;
            j++) {  //only check elements below the diagonal (col is lower than the row)
            double absErr =
                std::abs(FDHessian.coeff(i, j) - analyticHessian.coeff(i, j));
            double relError = 2 * absErr /
                              (eps + std::abs(FDHessian.coeff(i, j)) +
                               std::abs(analyticHessian.coeff(i, j)));
            if (relError > maxAcceptableError && absErr > maxAcceptableError) {
                Logger::print(Logger::RED,
                              "Mismatch element %d,%d: Analytic val: %lf, FD "
                              "val: %lf. Error: %lf(%lf%%)\n",
                              i, j, analyticHessian.coeff(i, j),
                              FDHessian.coeff(i, j), absErr, relError * 100);
            }
        }
    }
}

/*
void ObjectiveFunction::testSparseHessianPSD(const Vector &p) {
    / *
            if (svdAnalysis) {
                    SVDFactorization svd;
                    SparseMatrix U, S, V;
                    svd.computeSVD(*function->getHessianAt(pi), U, S, V);
                    for (int i = 0; i<min(S.getNumCols(), S.getNumRows()); i++)
                            Logger::printStatic("%d th singluer value of hessian
       : %e \n", i + 1, S.getElementAt(i, i));
            }
    * /

    SparseMatrix H(p.size(), p.size());
    / *
    DynamicArray<MTriplet> hessianEntries;
    hessianEntries.clear();
    addFiniteDifferenceHessianEntriesTo(p, hessianEntries);
    H.setFromTriplets(hessianEntries.begin(), hessianEntries.end());
* /
    computeHessian(p, H);

    Logger::print("Objective Function: testing hessian PDness...\n");

    Matrix Hd(H);
    Eigen::SelfAdjointEigenSolver<Matrix> es(Hd);
    Eigen::VectorXd D = es.eigenvalues();
    if (D.minCoeff() < 0)
        Logger::print(Logger::RED,
                      "Hessian is not PSD with min eigenvalues = %lf\n",
                      D.minCoeff());
    else
        Logger::print("Hessian is PSD with min eigenvalues = %lf\n",
                      D.minCoeff());
}
*/

//void ObjectiveFunction::computeHessian(const Vector &x, Matrix &hessian) const {
//    hessian.resize((int)x.size(), (int)x.size());
//    hessian.setZero();
//    addFiniteDifferenceHessianEntriesTo(x,hessian);
//}

//void ObjectiveFunction::addFiniteDifferenceHessianEntriesTo( const Vector &x, Matrix& hessian) const {
//    const double h = 1e-7;
//    for (int i = 0; i < x.size(); ++i) {
//        Vector dx(x.size());
//        dx.setZero();
//        dx[i] = h;

//        Vector gp(x.size()), gm(x.size());
//        computeGradient(x + dx, gp);
//        computeGradient(x - dx, gm);

//        Vector hess = (gp - gm) / (2. * h);
// only copy lower diagonal entries to the estimated hessian, as it is
// supposed to be sparse...
//        for (int j = i; j < x.size(); ++j) {
//            hessian(j, i) += hess[j];
//        }
//    }
//}

