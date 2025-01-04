#pragma once

#include <control/kinematics/EndEffector.h>
#include <control/GeneralizedCoordinatesRepresentation.h>

/**
 * A kinematic motion plan consists of a sequence of n (n >= 1) robot poses. Each pose
 * might have a preference for a particular target configuration (i.e. regularizer) which
 * can be stronger or weaker.
 */
class TargetPose {
public:
    int sIndex = 0;
    Vector target;
    Vector jw;

    TargetPose() {}

    TargetPose(const Vector &target, int sIndex, double w = 1.0) {
        this->sIndex = sIndex;
        this->target = target;
        jw.resize(target.size());
        jw.setOnes(); jw *= w;
    }
};

/**
 * This is the motion plan used in kinematic trajectory optimization tasks.
 */
class KinematicMotionPlan {

public:
    KinematicMotionPlan() {
        robot = nullptr;
    }

    KinematicMotionPlan(Robot* robot, int nSamples);

    ~KinematicMotionPlan(void) {}

    void getRobotPose(int index, RobotPose &rs);
    RobotPose getRobotPose(int index);
    void setRobotPoseFromMotionSample(int index);

    void writeOptimizationParametersToList(Vector &q);
    void setOptimizationParametersFromList(const Vector &q);

    inline int getStateDim() {
        return (freezeRootDOFs == false) ? robot->getJointCount() + 6
                                         : robot->getJointCount();
    }

    inline int getNumberOfMotionSamples() {
        return (int)robotPoses.size();
    }

    inline int getNumberOfOptimizationParameters(){
        return getStateDim() * getNumberOfMotionSamples();
    }

    void addEndEffector(const EndEffector& ee) {
        endEffectors.push_back(ee);
        endEffectors.back().qIdxOfParent = 6 + robot->idx(robot->getParentOf(ee.rigidBody));
    }

    void addTargetPose(const RobotPose& rp, int sIndex, double w = 1.0) {
        gcRR.setGeneralizedCoordinateValuesFromAFPose(rp);
        targetPoses.push_back(TargetPose(gcRR.getQ(), sIndex, w));
    }

    void setTargetPose(int tpIdx, const RobotPose& rp) {
        gcRR.setGeneralizedCoordinateValuesFromAFPose(rp);
        targetPoses[tpIdx].target = gcRR.getQ();
    }

    void setRobotPoseForMotionSample(const RobotPose& rp, int sIndex) {
        gcRR.setGeneralizedCoordinateValuesFromAFPose(rp);
        robotPoses[sIndex] = gcRR.getQ();
    }

public:
    //store a generalized coordinates representation, used during optimization
    GeneralizedCoordinatesRepresentation gcRR;

public:
    //as control objectives, we keep a list of end effector targets
    std::vector<EndEffector> endEffectors;

public:
    //the robot whose motion is represented by the kinematic motion plan
    Robot* robot;

    //the entire kinematic motion plan is represented through a sequence of robot poses
    Array<Vector> robotPoses;

    //and a list of pose regularizers
    Array<TargetPose> targetPoses;

    bool freezeRootDOFs = false;
    //if we're using 2nd order terms in the hessian of the EE objectives,
    //then derivatives will be correct, but the hessian may be indefinite,
    //so convergence suffers. Turn off for an approximate Hessian that is
    //always SPD.
    bool use2ndOrderEETerms = true;
};
