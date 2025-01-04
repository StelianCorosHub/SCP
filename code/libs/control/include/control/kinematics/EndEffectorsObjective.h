#pragma once

#include <control/kinematics/KinematicMotionPlan.h>

class EndEffectorsObjective {
    friend class IKUIHelper;
public:
    EndEffectorsObjective() { }
    EndEffectorsObjective(KinematicMotionPlan *mp, double weight = 1.0);
    ~EndEffectorsObjective(void) {}

    double evaluate() const;
    void addGradientTo(Vector& grad) const;
    void addHessianEntriesTo(std::vector<SMTriplet> &hessianEntries) const;

private:
    KinematicMotionPlan *kmPlan = nullptr;
    double w = 1.0;
};
