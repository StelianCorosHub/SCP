#pragma once

#include <control/kinematics/KinematicMotionPlan.h>
#include <control/kinematics/IKObjectiveFunction.h>
#include <optimization/NewtonMinimizer.h>
#include <RBSim/ArticulatedFigure.h>

class IKSolver {
public:
    IKSolver(Robot* robot, bool freezeRootDOFs = false);

    ~IKSolver(void) {}

    void solve(int nSteps = 10, bool useCurrentPoseAsRegularizer = false);

public:
    KinematicMotionPlan kmp;

    IKObjectiveFunction ikObjective;
    NewtonMinimizer minimizer;
    bool checkDerivatives = false;
};


