#include <RBSim/RBJoint.h>
#include <RBSim/RBUtils.h>
#include <utils/logger.h>
#include <utils/mathUtils.h>
#include <utils/utils.h>

int tmpRBJ = 0;

/**
 * Default constructor
 */
RBJoint::RBJoint(void) {
    Logger::print("creating an RBJoint: %d\n", ++tmpRBJ);
}

/**
 * Default destructor
 */
RBJoint::~RBJoint(void) {
    Logger::print("deleting an RBJoint: %d\n", --tmpRBJ);
}

void RBJoint::fixJointConstraints(bool fixChildPos, bool fixChildOrientation, bool fixChildLinearVel, bool fixChildAngularVel) {
    assert(child && parent);

    // first fix the orientation of the child given the orientation of the parent and the joint angle
    if (fixChildOrientation && isHingeJoint())
        child->orientation = parent->orientation * Quaternion(getJointAngle(), rotAxis);

    // now worry about the joint positions
    if (fixChildPos) {
        // compute the vector rc from the child's joint position to the child's
        // center of mass (in rbEngine coordinates)
        V3D rc = child->getWorldCoordinates(V3D(cJPos, P3D(0, 0, 0)));
        // and the vector rp that represents the same quanity but for the parent
        V3D rp = parent->getWorldCoordinates(V3D(pJPos, P3D(0, 0, 0)));

        // the location of the child's CM is now: pCM - rp + rc
        child->position = parent->position + (rc - rp);
    }

    if (fixChildAngularVel && isHingeJoint()) {
        child->angularVelocity = parent->angularVelocity +
                parent->getWorldCoordinates(rotAxis) * getJointSpeed();
    }

    if (fixChildLinearVel) {
        // we want to get the relative velocity at the joint to be 0. This can
        // be accomplished in many different ways, but this approach only
        // changes the linear velocity of the child
        V3D pJPosVel = parent->getVelocityForLocalCoordsPoint(pJPos);
        V3D cJPosVel = child->getVelocityForLocalCoordsPoint(cJPos);
        child->velocity -= cJPosVel - pJPosVel;
    }
}

