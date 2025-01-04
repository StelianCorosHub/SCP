#include <control/kinematics/KinematicMotionPlan.h>

KinematicMotionPlan::KinematicMotionPlan(Robot* robot, int nSamples) {
    this->robot = robot;

    gcRR = GeneralizedCoordinatesRepresentation(robot);

    for (int i = 0; i < nSamples; i++)
        robotPoses.push_back(gcRR.getQ());
}

void KinematicMotionPlan::getRobotPose(int index, RobotPose &rp) {
    gcRR.setQ(robotPoses[index]);
    rp = gcRR.getAFPose();
}

RobotPose KinematicMotionPlan::getRobotPose(int index){
    gcRR.setQ(robotPoses[index]);
    return gcRR.getAFPose();
}

void KinematicMotionPlan::setRobotPoseFromMotionSample(int index) {
    robot->setState(getRobotPose(index));
}

void KinematicMotionPlan::writeOptimizationParametersToList(Vector &q) {
    int stateDim = getStateDim();
    q.resize(stateDim * getNumberOfMotionSamples());

    // stack the robot poses into one array of optimization parameters
    for (int j = 0; j < getNumberOfMotionSamples(); j++) {
        for (int i = 0; i < stateDim; i++)
            if (freezeRootDOFs)
                q[j * stateDim + i] = robotPoses[j][i + 6];
            else
                q[j * stateDim + i] = robotPoses[j][i];
    }
}

void KinematicMotionPlan::setOptimizationParametersFromList(const Vector &q) {
    int stateDim = getStateDim();

    // unstack the robot states from one big array into the states at different
    // points in the motion plan
    for (int j = 0; j < getNumberOfMotionSamples(); j++) {
        for (int i = 0; i < stateDim; i++)
            if (freezeRootDOFs)
                robotPoses[j][i + 6] = q[j * stateDim + i];
            else
                robotPoses[j][i] = q[j * stateDim + i];
    }
}
