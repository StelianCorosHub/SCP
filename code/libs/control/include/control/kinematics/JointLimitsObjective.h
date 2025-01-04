#pragma once

#include <control/kinematics/KinematicMotionPlan.h>
#include <optimization/ObjectiveFunction.h>
#include <optimization/SoftUnilateralConstraint.h>
#include <RBSim/RBJoint.h>

class JointLimitsObjective {
    friend class IKUIHelper;
public:

    JointLimitsObjective() { }

    JointLimitsObjective(KinematicMotionPlan *kmPlan, double weight = 1.0) {
        this->kmPlan = kmPlan;
        this->w = weight;
    }

    ~JointLimitsObjective() {}

    double evaluate() const;
    void addGradientTo(Vector& grad) const;
    void addHessianEntriesTo(std::vector<SMTriplet> &hessianEntries) const;

public:
    KinematicMotionPlan *kmPlan = nullptr;
    double w = 1.0;
};
