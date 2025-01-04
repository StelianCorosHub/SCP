#pragma once

#include <RBSim/RBUtils.h>
#include <RBSim/RBSRepo.h>
#include <RBSim/AFState.h>
#include <RBSim/RBSRenderer.h>

#include <utils/basics.h>

/**
 * Articulated figures are hierarchical arrangements of rigid bodies that are connected
 * to each other via hinge (i.e. 1DOF) joints. Each articulated figure has one root, and every
 * rigid body has an arbitrary number of children. The root of the articulated figure has no parent,
 * and all other rigid bodies have exactly one parent. No kinematic loops are allowed in the
 * definition of an articulated figure.
 */
class ArticulatedFigure : public RBSRepo {
private:
    // to enable traversal of the articulated figure's hierarchy of rigid bodies and joints,
    // we will store, for each rigid body, its parent joint (nullptr for the root)
    Array<pRBJ> rbParents;
    // and a list of its child joints
    Array<Array<pRBJ>> rbChildren;

    void setupHierarchy() {
        rbIndices.clear();
        jIndices.clear();
        rbParents.clear();
        rbChildren.clear();

        Array<pRBJ> tmpJoints = joints;

        root = RBs[0];

        joints.clear();
        RBs.clear();

        RBs.push_back(root);
        rbParents.push_back(pRBJ());
        rbChildren.push_back(Array<pRBJ>());

        uint rbIdx = 0;

        // we will sort the list of RBs and joints in a way that is convenient for downstream processing.
        // In particular, the joint connecting a child RB to its parent will have a smaller index than
        // all children joints of the child RB.

        while (rbIdx < RBs.size()) {
            auto rb = RBs[rbIdx];
            rbIndices[rb] = rbIdx;

            // look through the full list of joints - whichever one has rb as a parent, we'll need to link it properly
            for (auto j : tmpJoints) {
                //found a child joint of rb
                if (j->parent == rb) {
                    //so, check if it's a valid joint first
                    if (j->isHingeJoint() == false)
                        throwError("Articulated figures only support 1DOF hinge joints. Joint %s violates this assumption.\n", j->name.c_str());

                    //make sure the child RB does not yet belong to the list of RBs of the articulated figure - this avoids all kinematic loops
                    for (uint i = 0; i < RBs.size(); i++)
                        if (RBs[i].get() == j->child.get())
                            throwError("Kinematic loop avoiding in the structure of an articulated figure -- not allowed.");

                    //passed the test - add j to the global list of joints in the articulated figure
                    jIndices[j] = (int)joints.size();
                    joints.push_back(j);
                    //make sure current rb lists this joint as one of its children
                    rbChildren[rbIdx].push_back(j);

                    //add the child RB to the list and link properly
                    RBs.push_back(j->child);
                    rbParents.push_back(j);
                    rbChildren.push_back(Array<pRBJ>());
                }
            }
            rbIdx++;
        }

        //now, fix all constraints imposed by the joints
        fixArticulationConstraints();
        initialRootPos = root->position;
        initialRootRot = root->orientation;
    }

protected:
    ArticulatedFigure() {}

    P3D initialRootPos;
    Quaternion initialRootRot;

public:
    // root configuration
    pRB root = pRB();

public:
    /** the constructor */
    ArticulatedFigure(const char* fName, bool loadVisuals = true);

    /** the destructor */
    virtual ~ArticulatedFigure(void){}

    ArticulatedFigure(const ArticulatedFigure& r) {
        printf("Error: Articulated figure copy constructor is not implemented\n");
        exit(0);
    }

    ArticulatedFigure& operator=(const ArticulatedFigure& r) {
        printf("Error: Articulated figure copy operator is not implemented");
	    exit(0);
    }

    // sets up the sorted list of RBs and joints, as well as the links needed for indexing
    // and hierarchical traversal of the structure of the articulated figure
    virtual void finalize() override {
        setupHierarchy();
    }

    pRBJ getParentOf(const pRB& rb) const {
        return rbParents[idx(rb)];
    }

    double getTotalMass(){
        // compute the mass of the robot
        double mass = root->Mass();
        for (auto j : joints)
            mass += j->child->Mass();
        return mass;
    }

    /**
     * this method is used to read the reduced state of the robot from the file
     */
    void loadStateFromFile(const char *fName);

    /**
     * this method is used to write the reduced state of the robot to a file
     */
    void writeStateToFile(const char *fName);

    /**
     * returns current state
     */
    AFState getState() const;

    /**
     * sets the state of the robot using the input
     */
    void setState(const AFState& state);

    /**
     * returns a zero state (i.e. all joint angles and velocities set to 0)
     */
    AFState getZeroState();

    /**
     * returns default state (all velocities set to 0, all joint angles set to their default values)
     * */
    AFState getDefaultState();

    /**
     * sets the state of the robot using default joint configurations
     */
    void setDefaultState() {
        setState(getDefaultState());
    }

    // sets pose (i.e. all velocities are zero)
    void setPose(const AFPose& state);

    // returns current pose of the robot
    AFPose getPose() const ;

    void translateTo(const P3D& p) {
        AFState s(this);
        s.rootPos = p;
        setState(s);
    }

    void setRandomState(double rVal) {
        AFState rs = getDefaultState();
        for (uint i = 0; i<rs.q.size(); i++){
            rs.q[i] += (getRandomNumberIn01Range() - 0.5) * rVal;
            rs.qDot[i] += (getRandomNumberIn01Range() - 0.5) * rVal * 0.1;
        }

        rs.rootPos.x += (getRandomNumberIn01Range() - 0.5) * rVal;
        rs.rootPos.y += (getRandomNumberIn01Range() - 0.5) * rVal;
        rs.rootPos.z += (getRandomNumberIn01Range() - 0.5) * rVal;

        rs.rootVel.x() += (getRandomNumberIn01Range() - 0.5) * rVal * 0.1;
        rs.rootVel.y() += (getRandomNumberIn01Range() - 0.5) * rVal * 0.1;
        rs.rootVel.z() += (getRandomNumberIn01Range() - 0.5) * rVal * 0.1;

        rs.rootQ.x() += (getRandomNumberIn01Range() - 0.5) * rVal;
        rs.rootQ.y() += (getRandomNumberIn01Range() - 0.5) * rVal;
        rs.rootQ.z() += (getRandomNumberIn01Range() - 0.5) * rVal;
        rs.rootQ.w() += (getRandomNumberIn01Range() - 0.5) * rVal;
        rs.rootQ.normalize();

        rs.rootAngVel.x() += (getRandomNumberIn01Range() - 0.5) * rVal * 0.1;
        rs.rootAngVel.y() += (getRandomNumberIn01Range() - 0.5) * rVal * 0.1;
        rs.rootAngVel.z() += (getRandomNumberIn01Range() - 0.5) * rVal * 0.1;
        setState(rs);
    }

    /**
     * sets the state of the robot using all-zero joint configurations
     */
    void setZeroState() {
        setState(getZeroState());
    }

    /**
     * changes the state of the body parts of the articulated figure such that
     * all joint constraints are satisfied. Note that in general there are many
     * projections on the constraint manifold defined by the joints, and this
     * implementation chooses an arbitrary one.
     */
    void fixArticulationConstraints(){
        AFState s = getState();
        setState(s);
    }

    void drawSkeleton(pRB rb) {
        RBSRenderer::drawSkeleton(rb, rbParents[idx(rb)], rbChildren[idx(rb)]);
    }

    virtual void draw() override {
        for (uint i = 0; i < RBs.size(); i++){
            pRB rb = RBs[i];

            // Draw skeleton view first
            if (df.showSkeleton)
                drawSkeleton(rb);

            // and joint information
            if (df.showJoints && rbParents[idx(rb)].get() != nullptr)
                RBSRenderer::drawJoint(rbParents[idx(rb)]);

            // Then draw collsion spheres
            if (df.showCollisionPrimitives)
                RBSRenderer::drawCollisionPrimitives(rb);

            // and end effectors
            if (df.showEndEffectors)
                RBSRenderer::drawEndEffectors(rb);

        }

        //now we draw the items that are transparent
        for (uint i = 0; i < RBs.size(); i++){
            pRB rb = RBs[i];

            // Then draw meshes (because of blending)
            if (df.showMeshes && df.showMOI == false) {
                if (df.showSkeleton || df.showJoints || df.showCollisionPrimitives || df.showEndEffectors || df.drawTransparent)
                    RBSRenderer::draw3DModels(rb, 0.2);
                else
                    RBSRenderer::draw3DModels(rb);
            }

            // and now MOIs
            if (df.showMOI)
                RBSRenderer::drawMOI(rb);
        }
    }

    void setControlMode(MotorControlMode mcm) {
        for (auto tmpJ : joints)
            tmpJ->motor.controlMode = mcm;
    }

    void setMotorsKpKd(double kp, double kd) {
        for (auto tmpJ : joints) {
            tmpJ->motor.kp = kp;
            tmpJ->motor.kd = kd;
        }
    }

    void setMotorTargets(const AFState& s) {
        for (int i = 0; i < getJointCount(); i++){
            joints[i]->motor.targetPosition = s.q[i];
            joints[i]->motor.targetSpeed = s.qDot[i];
        }
    }

    void addMotorImGuiMenuItems() {
        if (ImGui::TreeNode("Robot Control Options")) {
            //we will set parameters for all joints
            pRBJ j = joints[0];
            double mkp = j->motor.kp;
            double mkd = j->motor.kd;

            Array<const char*> items{"Motors Off", "Position Control", "Velocity Control", "Torque Control"}; // defined somewhere
            int selectedIndex = j->motor.controlMode; // you need to store this state somewhere
            const char* current_item = items[selectedIndex];

            if (ImGui::BeginCombo("##combo", current_item)) {
                for (int i = 0; i < items.size(); ++i) {
                    const bool isSelected = (selectedIndex == i);
                    if (ImGui::Selectable(items[i], isSelected)) {
                        selectedIndex = i;
                        current_item = items[selectedIndex];
                        setControlMode((MotorControlMode)(MotorControlMode::OFF + selectedIndex));
                    }

                    // Set the initial focus when opening the combo
                    // (scrolling + keyboard navigation focus)
                    if (isSelected) {
                        ImGui::SetItemDefaultFocus();
                    }
                }
                ImGui::EndCombo();
            }

            ImGui::InputDouble("Motor kp:", &mkp);
            ImGui::InputDouble("Motor kd:", &mkd);

            setMotorsKpKd(mkp, mkd);

            if (ImGui::Button("Set Default Pose")) {
                for (auto tmpJ : joints)
                    tmpJ->motor.targetPosition = tmpJ->defaultJointAngle;
            }

            if (ImGui::Button("Set Zero Pose")) {
                for (auto tmpJ : joints)
                    tmpJ->motor.targetPosition = 0;
            }
            ImGui::TreePop();
        }

    }

};

typedef std::shared_ptr<ArticulatedFigure> pAF;
typedef ArticulatedFigure Robot;
typedef AFState RobotState;
typedef AFPose RobotPose;
typedef std::shared_ptr<Robot> pRobot;


