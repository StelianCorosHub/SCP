#include <RBSim/ArticulatedFigure.h>
#include <RBSim/RBSLoader.h>

ArticulatedFigure::ArticulatedFigure(const char* fName, bool loadVisuals) {
    // load robot from rbLoader
    RBSLoader rbsLoader(fName, loadVisuals);
    rbsLoader.populate(this);
}

void ArticulatedFigure::setState(const AFState& state) {
    root->position = state.rootPos;
    root->orientation = state.rootQ;
    root->velocity = state.rootVel;
    root->angularVelocity = state.rootAngVel;

    for (uint j = 0; j < joints.size(); j++) {
        joints[j]->child->orientation = joints[j]->parent->orientation * Quaternion(state.q[j], joints[j]->rotAxis);
        joints[j]->child->angularVelocity = joints[j]->parent->angularVelocity +
                joints[j]->parent->getWorldCoordinates(joints[j]->rotAxis) * state.qDot[j];

        // and now set the linear position and velocity
        joints[j]->fixJointConstraints(true, false, true, false);
    }
}

/**
 * returns current state
 */
AFState ArticulatedFigure::getState() const {
    AFState afs((int)joints.size());

    afs.rootPos = root->position;
    afs.rootQ = root->orientation;
    afs.rootVel = root->velocity;
    afs.rootAngVel = root->angularVelocity;

    for (uint j = 0; j < joints.size(); j++) {
        afs.q[j] = joints[j]->getJointAngle();
        afs.qDot[j] = joints[j]->getJointSpeed();
    }

    return afs;
}

void ArticulatedFigure::setPose(const AFPose& state) {
    root->position = state.rootPos;
    root->orientation = state.rootQ;
    root->velocity = V3D();
    root->angularVelocity = V3D();

    for (uint j = 0; j < joints.size(); j++) {
        joints[j]->child->orientation = joints[j]->parent->orientation * Quaternion(state.q[j], joints[j]->rotAxis);
        joints[j]->child->angularVelocity = V3D();

        // and now set the linear position and velocity
        joints[j]->fixJointConstraints(true, false, false, false);
    }
}

/**
 * returns current state
 */
AFPose ArticulatedFigure::getPose() const {
    AFPose afp((int)joints.size());

    afp.rootPos = root->position;
    afp.rootQ = root->orientation;


    for (uint j = 0; j < joints.size(); j++)
        afp.q[j] = joints[j]->getJointAngle();

    return afp;
}

/**
 * returns a zero state (i.e. all joint angles and velocities set to 0)
 */
AFState ArticulatedFigure::getZeroState() {
    AFState afs((int)joints.size());

    //note: by default, AFState is initialized to all zeros, so not much to do here
    afs.rootPos = initialRootPos;

    return afs;
}

/**
 * returns default state (all velocities set to 0, all joint angles set to their default values)
 * */
AFState ArticulatedFigure::getDefaultState() {
    AFState afs((int)joints.size());

    afs.rootPos = initialRootPos;
    afs.rootQ = initialRootRot;

    for (uint j = 0; j < joints.size(); j++) {
        afs.q[j] = joints[j]->defaultJointAngle;
    }

    return afs;
}

void ArticulatedFigure::loadStateFromFile(const char *fName) {
    AFState state = getState();
    state.readFromFile(fName);
    setState(state);
}

void ArticulatedFigure::writeStateToFile(const char *fName) {
    AFState state = getState();
    state.writeToFile(fName);
}
