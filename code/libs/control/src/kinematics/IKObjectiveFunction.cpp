#include <control/kinematics/IKObjectiveFunction.h>

IKObjectiveFunction::IKObjectiveFunction(KinematicMotionPlan *kmPlan) {
    this->kmPlan = kmPlan;
    eeObjective = EndEffectorsObjective(kmPlan, 1.0);
    jLimitObjective = JointLimitsObjective(kmPlan, 10.0);
    poseRegularizer = PoseRegularizer(kmPlan, 0.01);
}

IKObjectiveFunction::~IKObjectiveFunction() {}

// this function is called before each and every optimization step starts --
// can be used to precompute various quantities, set regularizer parameters, etc...
void IKObjectiveFunction::prepareForOptimizationStep(const Vector &x) {
    x0 = x;
}

// this should always return the current value of the objective function
double IKObjectiveFunction::evaluate(const Vector &x) {
    kmPlan->setOptimizationParametersFromList(x);

    return eeObjective.evaluate() + jLimitObjective.evaluate() + poseRegularizer.evaluate();
}

void IKObjectiveFunction::computeGradient(const Vector &x, Vector &grad) {
    kmPlan->setOptimizationParametersFromList(x);
    resize(grad, x.size());
    eeObjective.addGradientTo(grad);
    jLimitObjective.addGradientTo(grad);
    poseRegularizer.addGradientTo(grad);
}

void IKObjectiveFunction::computeHessian(const Vector &x, SparseMatrix &hessian) {
    kmPlan->setOptimizationParametersFromList(x);
    hessian.resize(x.size(), x.size());
    hessian.setZero();
    triplets.clear();

    eeObjective.addHessianEntriesTo(triplets);
    jLimitObjective.addHessianEntriesTo(triplets);
    poseRegularizer.addHessianEntriesTo(triplets);

    hessian.setFromTriplets(triplets.begin(), triplets.end());
}

