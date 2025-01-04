#pragma once

#include <RBSim/RB.h>
#include <utils/mathUtils.h>

#include <RBSim/RBSimConstraint.h>

#include <string>
#include <unordered_map>

enum MotorControlMode {
    OFF = 0,     //torque = 0
    POSITION_CONTROL, //torque ~ kp * (desPos - pos) + kd (0 - vel)
    VELOCITY_CONTROL, //torque ~ kd * (desVel - vel)
    TORQUE_CONTROL   //torque = desTorque
};

class JMotor {
public:
    // joint control mode
    MotorControlMode controlMode = MotorControlMode::OFF;

    // control signal parameters
    double targetTorque = 0;
    double targetPosition = 0;
    double targetSpeed = 0;

    // gains for position/velocity/impedence control mode
    double kp = 12000.0;
    double kd = 0.1;

    //speed and torque limits
    double maxSpeed = 100.0;    // rad/s
    double maxTorque = INFINITY;  // N.m

    // in position or velocity mode, this is the torque actually applied
    double appliedMotorTorque = 0;

    bool hasTorqueLimits() {
        return maxTorque != INFINITY;
    }

};

/**
 * This class is used to implements hinge/revolute 1DOF joints that allow
 * relative rotation between the parent and the child bodies about a given axis
 */
class RBJoint {
public:
    // the name of the joint
    std::string name;
    // parent rigid body
    pRB parent = pRB();
    // this is the location of the joint on the parent - expressed in the
    // parent's local coordinates
    P3D pJPos = P3D(0, 0, 0);
    // this is the child link
    pRB child = pRB();
    // this is the location of the joint on the child - expressed in the child's
    // local coordinates
    P3D cJPos = P3D(0, 0, 0);
    // local coordinates of the rotation axis. Since the child and parent only
    // rotate relative to each other about this joint, the local coordinates for
    // the rotation axis are the same in parent and child frame
    V3D rotAxis = V3D(0, 0, 0);

    //joint limits and restitution coefficient
    double minAngle = -INFINITY, maxAngle = INFINITY, jLimitRestitutionCoeff = 0.2;

    // default joint angle
    double defaultJointAngle = 0;

    //control parameters
    JMotor motor;

    // selected in GUI
    bool selected = false;

    //parameters that control the stiffness of the joint
    double jERP = 0.2;
    double jCFR = 0;

public:
    /** Default constructor */
    RBJoint(void);

    /** Default destructor */
    ~RBJoint(void);

    /**
     * Returns the world position of the joint
     */
    P3D getWorldPosition(){
        return (child->getWorldCoordinates(cJPos) + parent->getWorldCoordinates(pJPos)) / 2.0;
    }

    /**
     * computes the relative orientation between the parent and the child rigid
     * bodies
     */
    Quaternion getRelativeOrientation(){
        // qc = qp * qRel
        return (parent->orientation.inverse() * child->orientation).normalized();
    }

    /**
     * This method is used to get the angular velocities of the child relative to its parent,
     * expressed in the parent's local coordinate frame.
     */
    V3D getLocalCoordsRelativeAngularVelocity() {
        return parent->getLocalCoordinates(child->angularVelocity - parent->angularVelocity);
    }

    /**
     * returns the angle this joint makes with its parent. If it is a hinge joint, the angle is in the
     * range -PI to PI about the rotation axis of the joint. Otherwise, we get the signed angle.
     */
    inline double getJointAngle() {
        if (isHingeJoint())
            return getRelativeOrientation().getRotationAngle(rotAxis);
        return getRelativeOrientation().getRotationAngle();
    }

    /**
     * returns the rate of change of the angle this joint makes with its parent. If it is a hinge joint,
     * the angular speed is measured about the rotation axis of the joint, and it can be positive or negative.
     * Otherwise, we get the magnitude of the relative angular velocity
     */
    inline double getJointSpeed() {
        if (isHingeJoint())
            return getLocalCoordsRelativeAngularVelocity().dot(rotAxis);
        return getLocalCoordsRelativeAngularVelocity().norm();
    }

    //joints can either be 1DOF hinge joints, or 3DOF ballAndSocket-type joints.
    //if a rotation axis is specified for this joint, then it is a hinge joint, otherwise
    //it is a ballAndSocket joint
    bool isHingeJoint() {
        return rotAxis.isZero() == false;
    }

    bool hasJointLimits() {
        return isHingeJoint() && IS_FINITE(minAngle) && IS_FINITE(maxAngle) && minAngle < maxAngle;
    }

    void applyFFMotorTorques(std::unordered_map<pRB, V3D>& rbTorques) {
        motor.appliedMotorTorque = 0;
        if (motor.controlMode == TORQUE_CONTROL) {
            motor.appliedMotorTorque = motor.targetTorque;
            V3D rotAxisW = child->getWorldCoordinates(rotAxis);
            rbTorques[parent] += rotAxisW * motor.targetTorque;
            rbTorques[child]  -= rotAxisW * motor.targetTorque;
        }
    }

    /**
     * This method is used to fix the errors in the joints (i.e. project state
     * of child such that joint configuration is consistent). The state of the
     * parent does not change.
     */
    void fixJointConstraints(bool fixChildPos, bool fixChildOrientation,
                             bool fixChildLinearVel, bool fixChildAngularVel);

    void addJointAnchorConstraints(Array<RBSimConstraint>& constraints) {
        //we want the location of the joint as seen from the coordinate frame of the two
        //rigid bodies to coincide in space. This means that in the absence of errors,
        //we want the velocity of the joint locations to match.

        Matrix3x3 ra = getSkewSymmetricMatrix(parent->getWorldCoordinates(V3D(P3D(), pJPos)));
        Matrix3x3 rb = getSkewSymmetricMatrix(child->getWorldCoordinates(V3D(P3D(), cJPos)));
        Matrix3x3 I = Matrix3x3::Identity();
        V3D vRel = parent->getVelocityForLocalCoordsPoint(P3D() + pJPos) - child->getVelocityForLocalCoordsPoint(P3D() + cJPos);
        V3D err(child->getWorldCoordinates(P3D() + cJPos), parent->getWorldCoordinates(P3D() + pJPos));

        for (int d = 0; d < 3; d++){
            constraints.push_back(RBSimConstraint(parent, child));
            constraints.back().l_i = I.row(d); constraints.back().a_i = -ra.row(d);
            constraints.back().l_j = I.row(d) * -1; constraints.back().a_j = rb.row(d);
            constraints.back().c = err[d];

            constraints.back().erp = jERP;
            constraints.back().cfr = jCFR;
        }
    }

    void addJointAxisConstraints(Array<RBSimConstraint>& constraints) {
        //we want here the relative orientation (w_i - w_j) to be unconstrained about the
        //rotation axis r, but fully constrained along any axis orthogonal to r. That means
        //that, given two vectors p, q orthonormal to r, the dot product between them and the
        //relative angular velocity should be 0: p*w_i - p*w_j = 0 & q*w_i - q*w_j = 0
        //As for the current value of the constraint, we need to look at the smallest rotation
        //that aligns the joint axes. This smallest rotation is obtained by a cross product,
        //whose magnitude sin(theta) is approximately theta for small rotation angles.
        //The relative angular velocity due to the constraint forces needs to make the cross
        //product 0.
        if (!isHingeJoint())
            return;

        V3D orthogonalAxes[2];
        getOrthonormalVectors(child->getWorldCoordinates(rotAxis), orthogonalAxes[0], orthogonalAxes[1]);

        V3D crossProduct = parent->getWorldCoordinates(rotAxis).cross(child->getWorldCoordinates(rotAxis));

        for (int d = 0; d < 2; d++) {
            constraints.push_back(RBSimConstraint(parent, child));
            constraints.back().l_i = V3D(); constraints.back().a_i = orthogonalAxes[d];
            constraints.back().l_j = V3D(); constraints.back().a_j = -orthogonalAxes[d];
            constraints.back().c = -crossProduct.dot(orthogonalAxes[d]);

            constraints.back().erp = jERP;
            constraints.back().cfr = jCFR;
        }
    }

    void pushRotAxisConstraint(Array<RBSimConstraint>& constraints) {
        V3D rotAxisW = child->getWorldCoordinates(rotAxis);

        constraints.push_back(RBSimConstraint(parent, child));
        constraints.back().l_i = V3D(); constraints.back().a_i = rotAxisW;
        constraints.back().l_j = V3D(); constraints.back().a_j = -rotAxisW;
    }

    void addJointMotorConstraints(Array<RBSimConstraint>& constraints, double sim_dt) {
        motor.appliedMotorTorque = 0;

        if (!isHingeJoint())
            return;

        if (motor.controlMode == POSITION_CONTROL || motor.controlMode == VELOCITY_CONTROL) {
            pushRotAxisConstraint(constraints);

            constraints.back().cf_val = &motor.appliedMotorTorque;

            if (motor.hasTorqueLimits()) {
                constraints.back().cf_min = -motor.maxTorque;
                constraints.back().cf_max = motor.maxTorque;
            }

            if (motor.controlMode == POSITION_CONTROL) {
                double dAngle = motor.targetPosition - getJointAngle();
                clamp(dAngle, -sim_dt * motor.maxSpeed, sim_dt * motor.maxSpeed);

                constraints.back().c = dAngle;

                // ERP = h*kp / (h*kp + kd)
                constraints.back().erp = sim_dt * motor.kp / (sim_dt * motor.kp + motor.kd);
                // CFM = 1 / (h*kp + kd)
                constraints.back().cfr = 1.0 / (sim_dt * motor.kp + motor.kd);
            } else {
//            Logger::consolePrint("%lf", getJointSpeed());
                constraints.back().cDotFF = -motor.targetSpeed;
                clamp(constraints.back().cDotFF, -motor.maxSpeed, motor.maxSpeed);
                // CFM = 1 / (h*kp + kd)
                constraints.back().cfr = 1.0 / (sim_dt * motor.kp + motor.kd);
            }
        }
    }

    //NOTE: motors and joint limits are somewhat redundant or fighting against each other.
    //      With regularizers and such, this does not seem to be a problem, but still, if numerical issues
    //      occur, we might have to make sure only one is active at a time, depending on the state of the joint
    void addJointLimitsConstraints(Array<RBSimConstraint>& constraints, double sim_dt) {
        if (hasJointLimits() == false)
            return;

        //check out the joint limits...
        double jAngle = getJointAngle(); //child q - parent q
        double jAngleDot = getJointSpeed();
        if (jAngle <= minAngle) {
            pushRotAxisConstraint(constraints);
            constraints.back().c = -(jAngle - minAngle);

            if (jAngleDot < 0)
                constraints.back().cDotFF = jAngleDot * jLimitRestitutionCoeff;

            constraints.back().erp = 0.2;
            constraints.back().cfr = 1e-5;
            //this constraint is only allowed to push away from the joint limit, not pull towards it...
            constraints.back().cf_max = 0;
        }

        if (jAngle >= maxAngle) {
            pushRotAxisConstraint(constraints);
            constraints.back().c = (maxAngle - jAngle);

            if (jAngleDot > 0)
                constraints.back().cDotFF = jAngleDot * jLimitRestitutionCoeff;

            constraints.back().erp = 0.2;
            constraints.back().cfr = 1e-5;
            //this constraint is only allowed to push away from the joint limit
            constraints.back().cf_min = 0;
        }
    }

    void addUnboundedConstraintsToList(Array<RBSimConstraint>& constraints, double sim_dt) {
        addJointAnchorConstraints(constraints);

        addJointAxisConstraints(constraints);

        if (motor.hasTorqueLimits() == false)
            addJointMotorConstraints(constraints, sim_dt);
    }

    void addBoundedConstraintsToList(Array<RBSimConstraint>& constraints, double sim_dt) {
        if (motor.hasTorqueLimits())
            addJointMotorConstraints(constraints, sim_dt);

        addJointLimitsConstraints(constraints, sim_dt);
    }

    void addMotorImGuiMenuItems() {
        if (ImGui::TreeNode("JMotor Control Options")) {

            Array<const char*> items{"Motors Off", "Position Control", "Velocity Control", "Torque Control"}; // defined somewhere
            int selectedIndex = motor.controlMode; // you need to store this state somewhere
            const char* current_item = items[selectedIndex];

            if (ImGui::BeginCombo("##combo", current_item)) {
                for (int i = 0; i < items.size(); ++i) {
                    const bool isSelected = (selectedIndex == i);
                    if (ImGui::Selectable(items[i], isSelected)) {
                        selectedIndex = i;
                        current_item = items[selectedIndex];
                            motor.controlMode = (MotorControlMode)(MotorControlMode::OFF + selectedIndex);
                    }

                    // Set the initial focus when opening the combo
                    // (scrolling + keyboard navigation focus)
                    if (isSelected) {
                        ImGui::SetItemDefaultFocus();
                    }
                }
                ImGui::EndCombo();
            }

            ImGui::InputDouble("Motor kp:", &motor.kp);
            ImGui::InputDouble("Motor kd:", &motor.kd);

            ImGui::InputDouble("target position:", &motor.targetPosition);
            ImGui::InputDouble("target velocity:", &motor.targetSpeed);

            ImGui::TreePop();
        }

    }

};

typedef std::shared_ptr<RBJoint> pRBJ;

inline pRBJ createJoint(const std::string name, const pRB& rbP, const pRB& rbC, const P3D& pjPos, const P3D& cjPos){
    pRBJ j = pRBJ(new RBJoint());
    j->parent = rbP;
    j->child = rbC;
    j->name = name;
    j->pJPos = pjPos;
    j->cJPos = cjPos;

    return j;
}
