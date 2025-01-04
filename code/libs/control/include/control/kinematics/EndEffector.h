#pragma once

#include <RBSim/RB.h>
#include <RBSim/ArticulatedFigure.h>
#include <utils/mathUtils.h>

class EndEffector {
public:
    // end effector belongs is defined on this rigid body...
    pRB rigidBody = nullptr;
    // this is the index of the parent joint of the RB in the coordinate frame of the parent
    int qIdxOfParent = -1;

    //and on this sample in the motion plan
    int sIndex = -1;

    P3D localCoordinates; //and its position/orientation, on this rigid body are
    P3D targetPosition;   //world coordinates target for the position of the end effector
    V3D localFrame[3];    //vectors in the EE-local coordinate frame which will have specified targets
    V3D targetFrame[3];   //target world-coordinate values for the vectors in the local frame above

    // if we want to control just the x,y or z component of the end effector,
    // use this mask...
    V3D positionMask = V3D(0.0, 0.0, 0.0);
    // same concept for
    V3D orientationMask = V3D(0.0, 0.0, 0.0);

    double posWeight = 1.0;
    double rotWeight = 1.0;

public:
    EndEffector() {}
    EndEffector(pRB rb, const P3D& localCoordinates = P3D(), int sIndex = 0) {
        this->rigidBody = rb;
        this->localCoordinates = localCoordinates;
        this->sIndex = sIndex;
    }

    ~EndEffector() {}

    void setLocalCoordinates(const P3D &ee) { localCoordinates = ee; }

    void setTargetPosition(const P3D &target) {
        targetPosition = target;
        positionMask = V3D(1.0, 1.0, 1.0);
    }

    void setTargetOrientation(const Quaternion &q) {
        for (int i = 0; i < 3; i++) {
            localFrame[i] = V3D();
            localFrame[i][i] = 1.0;
            setTargetDirectionForVector(localFrame[i], q * localFrame[i], i);
        }
    }

    void setTargetDirectionForVector(const V3D &vLocal, const V3D &vGlobal, int idx = 0) {
        localFrame[idx] = vLocal;
        targetFrame[idx] = vGlobal;
        orientationMask[idx] = 1.0;
    }

    void setTargetOrientationForLocalFrame(const Matrix3x3 &localQFrame, const Matrix3x3 &targetQFrame) {
        for (int i = 0; i < 3; i++)
            setTargetDirectionForVector(V3D(localQFrame.col(i)), V3D(targetQFrame.col(i)), i);
    }

    Matrix3x3 getTargetOrientation() const {
        Matrix3x3 targetMatrix;
        for (int i = 0; i < 3; i++) targetMatrix.col(i) = targetFrame[i];
        return targetMatrix;
    }

    bool constrainPosition() {
        return positionMask.norm() > 1e-10;
    }

    bool constrainOrientation() {
        return orientationMask.norm() > 1e-10;
    }
};
