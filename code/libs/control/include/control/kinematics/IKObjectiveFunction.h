#pragma once

#include <control/kinematics/KinematicMotionPlan.h>
#include <optimization/ObjectiveFunction.h>
#include <control/kinematics/EndEffectorsObjective.h>
#include <control/kinematics/JointLimitsObjective.h>
#include <control/kinematics/RegularizerObjective.h>

class IKObjectiveFunction : public ObjectiveFunction {
    friend class IKUIHelper;
public:
    IKObjectiveFunction() {}
    IKObjectiveFunction(KinematicMotionPlan *kmPlan);
    virtual ~IKObjectiveFunction(void);

    // this function is called before each and every optimization step starts --
    // can be used to precompute various quantities, set regularizer parameters, etc...
    virtual void prepareForOptimizationStep(const Vector &x) override;

    // this should always return the current value of the objective function
    virtual double evaluate(const Vector &x) override;

    virtual void computeGradient(const Vector &x, Vector &grad) override;

    virtual void computeHessian(const Vector &x, SparseMatrix &hessian) override;

private:
    KinematicMotionPlan *kmPlan = nullptr;
    //keep the starting guess for the optimization step -- updated after each optimization step
    Vector x0;
    //keep a list of triplets here to avoid needing to reallocate space with every optimization step
    Array<SMTriplet> triplets;

public:
    EndEffectorsObjective eeObjective;
    JointLimitsObjective jLimitObjective;
    PoseRegularizer poseRegularizer;
};
