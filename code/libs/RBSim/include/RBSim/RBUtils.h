#pragma once

#include <utils/geoms.h>
#include <utils/mathUtils.h>
#include <utils/utils.h>
#include <imgui.h>

#include <utils/utils.h>

enum RBS_KEYWORDS {
    RBS_RB = 0,
    RBS_END_RB,
    RBS_NAME,
    RBS_MASS,
    RBS_MOI,
    RBS_IS_STATIC,
    RBS_COLLISION_SPHERE,
    RBS_COLLISION_BOX,
    RBS_COLLISION_PLANE,
    RBS_MESH_NAME,
    RBS_MESH_COLOR,
    RBS_MESH_TRANSFORMATION,
    RBS_CHILD,
    RBS_PARENT,
    RBS_PPOS,
    RBS_CPOS,
    RBS_END_EFFECTOR,
    RBS_CONTROL_MODE,
    RBS_JOINT,
    RBS_JOINT_LIMITS,
    RBS_JOINT_MAX_SPEED,
    RBS_JOINT_MAX_TORQUE,
    RBS_JOINT_AXIS,
    RBS_DEFAULT_ANGLE,
    RBS_JOINT_END,
    RBS_FRICTION_COEFF,
    RBS_REST_COEFF,
    RBS_POSITION,
    RBS_VELOCITY,
    RBS_ORIENTATION,
    RBS_ANGULAR_VELOCITY,
    RBS_COMMENT,
    RBS_END_OF_FILE,
    RBS_UNKNOWN,
};

static Array<Keyword> RBSKeywords = {{"RB", RBS_RB},
                                     {"/End_RB", RBS_END_RB},
                                     {"name", RBS_NAME},
                                     {"mass", RBS_MASS},
                                     {"moi", RBS_MOI},
                                     {"collisionSphere", RBS_COLLISION_SPHERE},
                                     {"collisionPlane", RBS_COLLISION_PLANE},
                                     {"collisionBox", RBS_COLLISION_BOX},
                                     {"static", RBS_IS_STATIC},
                                     {"mesh", RBS_MESH_NAME},
                                     {"meshColor", RBS_MESH_COLOR},
                                     {"meshTransformation", RBS_MESH_TRANSFORMATION},
                                     {"child", RBS_CHILD},
                                     {"parent", RBS_PARENT},
                                     {"endEffector", RBS_END_EFFECTOR},
                                     {"jointPPos", RBS_PPOS},
                                     {"jointCPos", RBS_CPOS},
                                     {"RBJoint", RBS_JOINT},
                                     {"jointLimits", RBS_JOINT_LIMITS},
                                     {"jointMaxSpeed", RBS_JOINT_MAX_SPEED},
                                     {"jointMaxTorque", RBS_JOINT_MAX_TORQUE},
                                     {"jointAxis", RBS_JOINT_AXIS},
                                     {"defaultAngle", RBS_DEFAULT_ANGLE},
                                     {"controlMode", RBS_CONTROL_MODE},
                                     {"/End_Joint", RBS_JOINT_END},
                                     {"frictionCoefficient", RBS_FRICTION_COEFF},
                                     {"restitutionCoefficient", RBS_REST_COEFF},
                                     {"position", RBS_POSITION},
                                     {"velocity", RBS_VELOCITY},
                                     {"orientation", RBS_ORIENTATION},
                                     {"angularVelocity", RBS_ANGULAR_VELOCITY}
};

class DrawingFlagContainer{
public:
    // drawing flags
    bool showMeshes = true;
    bool showCollisionPrimitives = false;
    bool showEndEffectors = false;
    bool showJoints = false;
    bool showSkeleton = false;
    bool showMOI = false;
    bool drawTransparent = false;

    /**
     * Draw engine options to ImGui.
     */
    virtual void showDrawOptions() {
        ImGui::Checkbox("Draw Meshes", &showMeshes);
        ImGui::Checkbox("Draw Skeleton", &showSkeleton);
        ImGui::Checkbox("Draw Joints", &showJoints);
        ImGui::Checkbox("Draw Collision Primitives", &showCollisionPrimitives);
        ImGui::Checkbox("Draw MOI box", &showMOI);
        ImGui::Checkbox("Draw End Effectors", &showEndEffectors);
        ImGui::Checkbox("Draw Transparent", &drawTransparent);
    }
};

/**
 * returns the angular velocity that explains how we got from qStart to qEnd in
 * dt time (qEnd = qDueToAngVelOverDT * qStart)
 */
inline V3D estimateAngularVelocity(const Quaternion &qStart, const Quaternion &qEnd, double dt) {
    // qEnd = rot(w_p, dt) * qStart
    Quaternion qRot = qEnd * qStart.inverse();

    V3D rotAxis(qRot.vec());
    if (rotAxis.norm() < 1e-10)
        return V3D(0, 0, 0);
    rotAxis.normalize();
    double rotAngle = (qRot.getRotationAngle());

    // qRot is the result of rotating with this fixed angular velocity for
    // some time dt...
    return rotAxis * rotAngle / dt;
}

inline Quaternion integrateOrientationForwardInTime(const Quaternion& qStart, const V3D& angularVelocity, double dt) {
    //qEnd = rot(w_p, dt) * qStart
    return Quaternion(dt * angularVelocity.norm(), angularVelocity.unit()) * qStart;
}
