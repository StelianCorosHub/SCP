#pragma once

#include <control/kinematics/KinematicMotionPlan.h>
#include <optimization/ObjectiveFunction.h>

/**
 * Regularizer: Don't go too far away from nominal pose
 */
class PoseRegularizer {
    friend class IKUIHelper;
public:
    PoseRegularizer() { }
    PoseRegularizer(KinematicMotionPlan *mp, double weight = 1.0);
    ~PoseRegularizer(void) {}

    double evaluate() const;
    void addGradientTo(Vector& grad) const;
    void addHessianEntriesTo(std::vector<SMTriplet> &hessianEntries) const;

private:
    KinematicMotionPlan *kmPlan;
public:
    double w = 1.0;
};

/*
/ **
 * Regularizer: penalize large joint angle accelerations
 * /
class SmoothMotionRegularizer {
public:
    SmoothMotionRegularizer(KinematicMotionPlan *mp, int startQIndex,
                            int endQIndex, std::string description,
                            double weight);
    ~SmoothMotionRegularizer(void) {}

    double evaluate() const;
    void addGradientTo(Vector& grad) const;
    void addHessianEntriesTo(std::vector<SMTriplet> &hessianEntries) const;

private:
    KinematicMotionPlan *kmPlan;
    int startQIndex = 0, endQIndex = 0;
};

/ **
 * Regularizer: have start and end states be specified, and 0 velocity
 * at the start and end"
 * /
class TargetStartAndEndPoseRegularizer {
public:
    TargetStatesObjective(KinematicMotionPlan *mp, std::string description,
                          double weight);
    ~TargetStatesObjective(void) {}

    double evaluate() const;
    void addGradientTo(Vector& grad) const;
    void addHessianEntriesTo(std::vector<SMTriplet> &hessianEntries) const;

private:
    KinematicMotionPlan *kmPlan = nullptr;
};
*/
