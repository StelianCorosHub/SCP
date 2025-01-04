#include <optimization/SoftUnilateralConstraint.h>


SoftLowerLimitConstraint::SoftLowerLimitConstraint(double l, double stiffness,
                                                   double epsilon) {
    this->limit = l;
    setPolynomialCoeffs(stiffness, epsilon);
}

double SoftLowerLimitConstraint::getEpsilon(){
    return epsilon;
}

double SoftLowerLimitConstraint::getStiffness(){
    return a1;
}

//based on the value of the stiffness and eps parameters, update the coeffs
//of the piece-wise cubic polynomial that we're working with here
void SoftLowerLimitConstraint::setPolynomialCoeffs(double stiffness, double eps){
    this->epsilon = eps;
    a1 = stiffness;
    b1 = -0.5 * a1 * epsilon;
    c1 = -1.0 / 3 * (-b1 - a1 * epsilon) * epsilon -
         1.0 / 2 * a1 * epsilon * epsilon - b1 * epsilon;

    a2 = (-b1 - a1 * epsilon) / (epsilon * epsilon);
    b2 = a1;
    c2 = b1;
    d2 = c1;
}

SoftLowerLimitConstraint::~SoftLowerLimitConstraint() {}

void SoftLowerLimitConstraint::setLimit(double l) { limit = l; }

double SoftLowerLimitConstraint::getLimit(){
    return limit;
}

// returns 1/2 C'C, where C is the current set of equality constraint values
double SoftLowerLimitConstraint::evaluate(double x) {
    x = x - limit;
    if (x < 0) return 0.5 * a1 * x * x + b1 * x + c1;
    if (x < epsilon)
        return 1.0 / 3 * a2 * x * x * x + 0.5 * b2 * x * x + c2 * x + d2;
    return 0;
}

double SoftLowerLimitConstraint::computeDerivative(double x) {
    x = x - limit;
    if (x < 0) return a1 * x + b1;
    if (x < epsilon) return a2 * x * x + b2 * x + c2;
    return 0;
}

double SoftLowerLimitConstraint::computeSecondDerivative(double x) {
    x = x - limit;
    if (x < 0) return a1;
    if (x < epsilon) return 2 * a2 * x + b2;
    return 0;
}

SoftUpperLimitConstraint::SoftUpperLimitConstraint(double u, double stiffness,
                                                   double epsilon) {
    this->limit = u;
    setPolynomialCoeffs(stiffness, epsilon);
}

SoftUpperLimitConstraint::~SoftUpperLimitConstraint() {}

void SoftUpperLimitConstraint::setLimit(double l) { limit = l; }

//based on the value of the stiffness and eps parameters, update the coeffs
//of the piece-wise cubic polynomial that we're working with here
void SoftUpperLimitConstraint::setPolynomialCoeffs(double stiffness, double eps){
    this->epsilon = eps;
    a1 = stiffness;

    b1 = 0.5 * a1 * epsilon;
    c1 = 1. / 6. * a1 * epsilon * epsilon;
    a2 = 1. / (2. * epsilon) * a1;
    b2 = a1;
    c2 = 0.5 * a1 * epsilon;
    d2 = 1. / 6. * a1 * epsilon * epsilon;
}

double SoftUpperLimitConstraint::getEpsilon(){
	return epsilon;
}

double SoftUpperLimitConstraint::getStiffness(){
	return a1;
}

double SoftUpperLimitConstraint::getLimit(){
    return limit;
}

// returns 1/2 C'C, where C is the current set of equality constraint values
double SoftUpperLimitConstraint::evaluate(double x) {
    x = x - limit;
    if (x > 0) return 0.5 * a1 * x * x + b1 * x + c1;
    if (x > -epsilon)
        return 1.0 / 3 * a2 * x * x * x + 0.5 * b2 * x * x + c2 * x + d2;
    return 0;
}

double SoftUpperLimitConstraint::computeDerivative(double x) {
    x = x - limit;
    if (x > 0) return a1 * x + b1;
    if (x > -epsilon) return a2 * x * x + b2 * x + c2;
    return 0;
}

double SoftUpperLimitConstraint::computeSecondDerivative(double x) {
    x = x - limit;
    if (x > 0) return a1;
    if (x > -epsilon) return 2 * a2 * x + b2;
    return 0;
}

SoftBoundConstraint::SoftBoundConstraint(double lLimit, double uLimit,
                                         double stiffness, double relEps)
    : l(lLimit, stiffness, fabs(uLimit - lLimit) * relEps),
      u(uLimit, stiffness, fabs(uLimit - lLimit) * relEps) {}

SoftBoundConstraint::~SoftBoundConstraint() {}

double SoftBoundConstraint::getEPS(){
    return l.getEpsilon() / fabs(u.getLimit() - l.getLimit());
}

double SoftBoundConstraint::getStiffness(){
	return l.getStiffness();
}

void SoftBoundConstraint::setStiffnessAndRelEpsilon(double stiffness, double eps) {
    l.setPolynomialCoeffs(stiffness, fabs(u.getLimit() - l.getLimit()) * eps);
    u.setPolynomialCoeffs(stiffness, fabs(u.getLimit() - l.getLimit()) * eps);
}

double SoftBoundConstraint::evaluate(double x) {
    return l.evaluate(x) + u.evaluate(x);
}

double SoftBoundConstraint::computeDerivative(double x) {
    return l.computeDerivative(x) + u.computeDerivative(x);
}

double SoftBoundConstraint::computeSecondDerivative(double x) {
    return l.computeSecondDerivative(x) + u.computeSecondDerivative(x);
}
