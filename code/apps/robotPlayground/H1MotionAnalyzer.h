#pragma once

#include <gui/application.h>

#include <utils/BVHReader.h>

#include "H1.h"
#include <control/kinematics/IKSolver.h>

#include "MocapBindings.h"

/**
 * Motion Analyzer for humanoid robots -> The assumption is that the BVH mocap clips follow the structrure from the 100styles dataset.
 **/
class H1MotionAnalyzer {
public:

    ~H1MotionAnalyzer() {
        delete ikSolver;
    }

    void initialize(pH1 h1, BVHClip* mocapClip) {
        this->h1 = h1;
        this->mocapClip = mocapClip;
        
        mocapClip->skeleton.scale = h1->getHeight() / mocapClip->skeleton.getUnscaledHeight();

        h1->setZeroState();
        MocapSkeletonPose defaultPose = mocapClip->getDefaultSkeletonPose();
        defaultPose.translateTo(h1->root->position);

        heightOffset = defaultPose.getBoneEndPos(R_TOES_IDX).y - h1->getLeftHeelPosition().y;

        initialMocapRootPosition = mocapClip->getSkeletonPoseAtFrame(0).getRootPosition();
        initialMocapRootPosition.y = 0;

        ikSolver = new IKSolver(h1.get(), true);
        ikSolver->minimizer.reg = 0.1;
        ikSolver->kmp.use2ndOrderEETerms = true;

        ikSolver->kmp.targetPoses[0].jw *= 0;
        ikSolver->kmp.targetPoses[0].jw[6 + h1->idx(h1->getJointByName("left_shoulder_yaw_joint"))] = 100;
        ikSolver->kmp.targetPoses[0].jw[6 + h1->idx(h1->getJointByName("left_elbow_pitch_joint"))] = 1;
        ikSolver->kmp.targetPoses[0].jw[6 + h1->idx(h1->getJointByName("right_shoulder_yaw_joint"))] = 100;
        ikSolver->kmp.targetPoses[0].jw[6 + h1->idx(h1->getJointByName("right_elbow_pitch_joint"))] = 1;

        ikSolver->kmp.addEndEffector(EndEffector(h1->getRBByName("left_hip_roll_link"), P3D()));
        ikSolver->kmp.addEndEffector(EndEffector(h1->getRBByName("left_knee_link"), P3D()));
        ikSolver->kmp.addEndEffector(EndEffector(h1->getRBByName("left_ankle_roll_link"), P3D()));

        ikSolver->kmp.addEndEffector(EndEffector(h1->getRBByName("right_hip_roll_link"), P3D()));
        ikSolver->kmp.addEndEffector(EndEffector(h1->getRBByName("right_knee_link"), P3D()));
        ikSolver->kmp.addEndEffector(EndEffector(h1->getRBByName("right_ankle_roll_link"), P3D()));

        ikSolver->kmp.addEndEffector(EndEffector(h1->getRBByName("torso_link"), P3D()));

        ikSolver->kmp.addEndEffector(EndEffector(h1->getRBByName("left_shoulder_yaw_link"), P3D()));
        ikSolver->kmp.addEndEffector(EndEffector(h1->getRBByName("left_elbow_pitch_link"), P3D()));

        ikSolver->kmp.addEndEffector(EndEffector(h1->getRBByName("right_shoulder_yaw_link"), P3D()));
        ikSolver->kmp.addEndEffector(EndEffector(h1->getRBByName("right_elbow_pitch_link"), P3D()));
    }

    MocapSkeletonPose getMocapPose(int fIdx) {
        MocapSkeletonPose msp = mocapClip->getDefaultSkeletonPose();
        msp = mocapClip->getSkeletonPoseAtFrame(fIdx);
        msp.translateBy(V3D(initialMocapRootPosition, P3D(0, 0, 0)));
        return msp;
    }

    MocapSkeletonPose getDefaultMocapPose() {
        MocapSkeletonPose msp = mocapClip->getDefaultSkeletonPose();
        msp.translateBy(V3D(initialMocapRootPosition, P3D(0, 0, 0)) + V3D(0, 0.91, 0));
        return msp;
    }

    inline Matrix3x3 getRot(const V3D& x, const V3D& y) {
        Matrix3x3 res;
        res.col(0) = x;
        res.col(1) = y;
        res.col(2) = x.cross(y);
        return res;
    }

    V3D getLegPlaneNormal(const Quaternion& rootQ, const V3D& hipToKneeV, const V3D& kneeToAnkleV, const Quaternion& footQ){
        V3D legPlaneNormal = hipToKneeV.cross(kneeToAnkleV).unit();
        double legAngle = hipToKneeV.angleWith(kneeToAnkleV);

        if (legPlaneNormal.dot(rootQ * H1_RIGHT) < 0)
            legAngle *= -1;

        //when the leg is nearly straight, we'll be using the foot longitudinal axis to determine the normal
        //of the plane that the limb lies within
        double t = mapTo01Range(legAngle, RAD(5), RAD(15));
        legPlaneNormal = legPlaneNormal * t + -hipToKneeV.cross(footQ * H1_FORWARD).unit() * (1-t);
        legPlaneNormal.normalize();

        return legPlaneNormal;
    }

    double getAnkleFEAngle(const Matrix3x3& R_lowerLeg, const V3D& ankleToToesV, double defaultShinFootAngle) {
        //bring the ankle-toes vector to the coordinate frame of the lower leg
        V3D localAnkleToToesV = V3D(R_lowerLeg.transpose() * ankleToToesV);

        //then project out the component that doesn't look like a flexion/extension motion
        localAnkleToToesV.setComponentAlong(H1_RIGHT, 0);

        //now compute the angle between the shin and this vector, minus the default angle
        return -defaultShinFootAngle + V3D(0,-1,0).signedAngleWith(localAnkleToToesV, H1_RIGHT);
    }

    double getToesAngle(const Matrix3x3& R_lowerLeg, const V3D& ankleToToesV, const V3D& toesV, double defaultFootToesAngle){
        return -defaultFootToesAngle + ankleToToesV.signedAngleWith(toesV, R_lowerLeg * H1_RIGHT);
    }

    void setH1PoseFromMocapSkeletonPose(const MocapSkeletonPose& msp) {
        h1->setDefaultState();
        H1State s(h1);
        s.rootPos = msp.getSFPos(0) + UP * heightOffset;
        s.rootQ = msp.getRootOrientation();
        h1->setState(s);

        V3D lHipToKneeV = msp.getBoneVector(L_UPPERLEG_IDX);
        //Logger::consolePrint("l upper leg: %lf %lf %lf\n", lHipToKneeV.x(), lHipToKneeV.y(), lHipToKneeV.z());
        ikSolver->kmp.endEffectors[0].setTargetDirectionForVector(V3D(0, -0.435985, -0.000061), lHipToKneeV);

        V3D lKneeToAnkleV = msp.getBoneVector(L_LOWERLEG_IDX);
        //Logger::consolePrint("l lower leg: %lf %lf %lf\n", lKneeToAnkleV.x(), lKneeToAnkleV.y(), lKneeToAnkleV.z());
        ikSolver->kmp.endEffectors[1].setTargetDirectionForVector(V3D(0, -0.430705, -0.000896), lKneeToAnkleV);

        V3D lAnkleToToesV = msp.getBoneVector(L_FOOT_IDX);
        //Logger::consolePrint("l foot: %lf %lf %lf\n", lAnkleToToesV.x(), lAnkleToToesV.y(), lAnkleToToesV.z());
        ikSolver->kmp.endEffectors[2].setTargetDirectionForVector(V3D(0, -0.096676, 0.190879), lAnkleToToesV);

        V3D rHipToKneeV = msp.getBoneVector(R_UPPERLEG_IDX);
        //Logger::consolePrint("r upper leg %lf %lf %lf\n", rHipToKneeV.x(), rHipToKneeV.y(), rHipToKneeV.z());
        ikSolver->kmp.endEffectors[3].setTargetDirectionForVector(V3D(0, -0.435985, -0.000061), rHipToKneeV);

        V3D rKneeToAnkleV = msp.getBoneVector(R_LOWERLEG_IDX);
        //Logger::consolePrint("r lower leg %lf %lf %lf\n", rKneeToAnkleV.x(), rKneeToAnkleV.y(), rKneeToAnkleV.z());
        ikSolver->kmp.endEffectors[4].setTargetDirectionForVector(V3D(0, -0.430705, -0.000896), rKneeToAnkleV);

        V3D rAnkleToToesV = msp.getBoneVector(R_FOOT_IDX);
        //Logger::consolePrint("r foot %lf %lf %lf\n", rAnkleToToesV.x(), rAnkleToToesV.y(), rAnkleToToesV.z());
        ikSolver->kmp.endEffectors[5].setTargetDirectionForVector(V3D(0, -0.096676, 0.190879), rAnkleToToesV);

        V3D torsoV = msp.getBoneVector(TORSO_IDX);
        //Logger::consolePrint("torso up %lf %lf %lf\n", torsoV.x(), torsoV.y(), torsoV.z());
        ikSolver->kmp.endEffectors[6].setTargetDirectionForVector(V3D(0, 0.132521, 0), torsoV, 0);

        V3D shoulderV = V3D(msp.getBoneStartPos(R_UPPER_ARM_IDX), msp.getBoneStartPos(L_UPPER_ARM_IDX));
        //Logger::consolePrint("torso s %lf %lf %lf\n", shoulderV.x(), shoulderV.y(), shoulderV.z());
        ikSolver->kmp.endEffectors[6].setTargetDirectionForVector(V3D(0.378450, 0, 0), shoulderV, 1);

        V3D lUpperArmV = msp.getBoneVector(L_UPPER_ARM_IDX);
        //Logger::consolePrint("l upper arm %lf %lf %lf\n", lUpperArmV.x(), lUpperArmV.y(), lUpperArmV.z());
        ikSolver->kmp.endEffectors[7].setTargetDirectionForVector(V3D(0, -0.310483, 0), lUpperArmV);

        V3D lLowerArmV = msp.getBoneVector(L_LOWER_ARM_IDX);
        //Logger::consolePrint("l lower arm %lf %lf %lf\n", lLowerArmV.x(), lLowerArmV.y(), lLowerArmV.z());
        ikSolver->kmp.endEffectors[8].setTargetDirectionForVector(V3D(0, 0, 0.254514), lLowerArmV);

        V3D rUpperArmV = msp.getBoneVector(R_UPPER_ARM_IDX);
        //Logger::consolePrint("r upper arm %lf %lf %lf\n", rUpperArmV.x(), rUpperArmV.y(), rUpperArmV.z());
        ikSolver->kmp.endEffectors[9].setTargetDirectionForVector(V3D(0, -0.310483, 0), rUpperArmV);

        V3D rLowerArmV = msp.getBoneVector(R_LOWER_ARM_IDX);
        //Logger::consolePrint("r lower arm %lf %lf %lf\n", rLowerArmV.x(), rLowerArmV.y(), rLowerArmV.z());
        ikSolver->kmp.endEffectors[10].setTargetDirectionForVector(V3D(0, 0, 0.254514), rLowerArmV);

        ikSolver->solve(200);

        s = H1State(h1);

        V3D v = h1->getLeftFootOrientation() * H1_RIGHT;
        double angleErrBefore = v.y();
        V3D rollAxis = h1->getLeftFootOrientation() * H1_FORWARD;
        V3D vProj = v; vProj.y() = 0; vProj.setComponentAlong(rollAxis, 0); //vProj.normalize();
        s.q[h1->getLeftAnkleRollIndex()] += v.signedAngleWith(vProj, rollAxis);

//        Logger::consolePrint("vProj.cross(v).angleWith(rollAxis): %lf, vProj.y(): %lf\n", vProj.cross(v).angleWith(rollAxis), vProj.y());

        v = h1->getRightFootOrientation() * H1_RIGHT;
        angleErrBefore = v.y();
        rollAxis = h1->getRightFootOrientation() * H1_FORWARD;
        vProj = v; vProj.y() = 0; vProj.setComponentAlong(rollAxis, 0); //vProj.normalize();
        s.q[h1->getRightAnkleRollIndex()] += v.signedAngleWith(vProj, rollAxis);

        h1->setState(s);

//        v = h1->getLeftFootOrientation() * H1_RIGHT;
//        Logger::consolePrint("angle err before: %lf, after: %lf\n", angleErrBefore, v.y());
    }

public:
    pH1 h1;
    BVHClip* mocapClip;

    P3D initialMocapRootPosition;

    //this is to make sure that H1's feet are on the ground in stance
    double heightOffset = 0;

    IKSolver *ikSolver = nullptr;
};

