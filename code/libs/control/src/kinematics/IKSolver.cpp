#include <control/kinematics/IKSolver.h>
#include <utils/logger.h>

IKSolver::IKSolver(Robot* robot, bool freezeRootDOFs) {
    kmp = KinematicMotionPlan(robot, 1);
    kmp.freezeRootDOFs = freezeRootDOFs;
    //add one pose regularizer...
    kmp.addTargetPose(RobotPose(robot), 0);
    ikObjective = IKObjectiveFunction(&kmp);
    minimizer.reg = 100;
    minimizer.printOutput = false;
}

void IKSolver::solve(int nSteps, bool useCurrentPoseAsRegularizer) {
    // will start the optimization from the current pose of the robot
    kmp.setRobotPoseForMotionSample(RobotPose(kmp.robot), 0);

    // resets the pose used as a regularizer for the IK solver
    if (useCurrentPoseAsRegularizer)
        kmp.setTargetPose(0, RobotPose(kmp.robot));

    Vector x;
    kmp.writeOptimizationParametersToList(x);

    if (checkDerivatives) {
        ikObjective.testGradientWithFD(x);
        ikObjective.testSparseHessianWithFD(x);
    }

    minimizer.minimize(&ikObjective, x, nSteps);

    if (minimizer.printOutput)
        Logger::consolePrint("Final IK energy val: %lf\n",
                             ikObjective.evaluate(x));

    kmp.setRobotPoseFromMotionSample(0);
}
