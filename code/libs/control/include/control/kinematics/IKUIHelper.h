#pragma once

#include <gui/application.h>
#include <gui/camera.h>
#include <gui/light.h>
#include <gui/renderer.h>
#include <gui/shader.h>
#include <control/kinematics/IKSolver.h>
#include <RBSim/ArticulatedFigure.h>
#include <utils/logger.h>

/**
 * UI helper class for IK. Manages the setting of end effector targets, drawing debug information, ImGui menus, etc
 *
 * Left Click + drag: set and move target on selected robot body part
 * Middle Click (or left CTRL + left click): remove last end effector target
 * Right Click: bring up the translate/rotate widget that will set end effector
 * targets for positions/orientations...
 */

class IKUIHelper {
public:
    IKUIHelper(IKSolver* ikSolver){
        this->ikSolver = ikSolver;
    }

    IKUIHelper(){
        this->ikSolver = nullptr;
    }

    ~IKUIHelper(){}

    bool mouseDrag(MouseState mouseState, const Ray& mouseRay) {
        if (mouseState.lButtonPressed == true) {
            if (selectedRB != NULL) {
                mouseRay.getDistanceToPoint(selectedPoint, &targetForSelectedPoint);
                ikSolver->kmp.endEffectors.back().setTargetPosition(targetForSelectedPoint);

                return selectedRB != NULL;
            }
        }
        return false;
    }

    bool mouseButtonReleased(int button, int mods){
        if (button == GLFW_MOUSE_BUTTON_LEFT && selectedRB) {
            selectedRB->selected = false;
            selectedRB = nullptr;
            if (usePersistentEETargets == false)
                ikSolver->kmp.endEffectors.clear();
        }
        return true;
    }

    bool mouseButtonPressed(int button, int mods, const Ray& mouseRay, int nEEMax = 0) {
        if (button == GLFW_MOUSE_BUTTON_MIDDLE || (button == GLFW_MOUSE_BUTTON_LEFT && mods & GLFW_MOD_CONTROL)){
            if (ikSolver->kmp.endEffectors.size() > (uint)nEEMax)
                ikSolver->kmp.endEffectors.pop_back();
            return true;
        }

        if (button == GLFW_MOUSE_BUTTON_LEFT || button == GLFW_MOUSE_BUTTON_RIGHT) {
            if (selectedRB != NULL)
                selectedRB->selected = false;
            showTargetSelectionGizmo = false;

            selectedRB = ikSolver->kmp.robot->getFirstRBHitByRay(mouseRay, selectedPoint);

            if (selectedRB != NULL) {
                Logger::consolePrint("clicked on %s at %lf %lf %lf\n", selectedRB->name.c_str(),
                                     selectedRB->getLocalCoordinates(selectedPoint).x,
                                     selectedRB->getLocalCoordinates(selectedPoint).y,
                                     selectedRB->getLocalCoordinates(selectedPoint).z);

                selectedRB->selected = true;
                targetForSelectedPoint = selectedPoint;
                targetOrientationForSelectedRB = selectedRB->orientation;
                // if we're using the gizmos, end effector targets get set when
                // we know if we want rotations or translations of the end
                // effector...
                if (button == GLFW_MOUSE_BUTTON_RIGHT) {
                    showTargetSelectionGizmo = true;
                    addedEE = false;
                    oldT = targetForSelectedPoint;
                    oldQ = targetOrientationForSelectedRB;
                } else {
                    ikSolver->kmp.addEndEffector(EndEffector(selectedRB, selectedRB->getLocalCoordinates(
                            targetForSelectedPoint)));
                    ikSolver->kmp.endEffectors.back().setTargetPosition(
                        targetForSelectedPoint);
                }
            }
            return true;
        }
        return false;
    }


    // objects drawn here will not have shadows cast on them
    void drawDebugInfo(){
        if (selectedRB)
            drawSphere(targetForSelectedPoint, 0.05f);

        for (uint i = 0; i < ikSolver->kmp.endEffectors.size(); i++) {
            if (ikSolver->kmp.endEffectors[i].constrainPosition()) {
                P3D p1 =
                    ikSolver->kmp.endEffectors[i]
                        .rigidBody->getWorldCoordinates(
                            ikSolver->kmp.endEffectors[i].localCoordinates);
                P3D p2 = ikSolver->kmp.endEffectors[i].targetPosition;
                drawCapsule(p1, p2, 0.025);
            }
            if (ikSolver->kmp.endEffectors[i].constrainOrientation()) {
                P3D p =
                    ikSolver->kmp.endEffectors[i]
                        .rigidBody->getWorldCoordinates(
                            ikSolver->kmp.endEffectors[i].localCoordinates);
                for (int j = 0; j < 3; j++) {
                    V3D v1 = 0.2 * ikSolver->kmp.endEffectors[i]
                                       .rigidBody->getWorldCoordinates(
                                           ikSolver->kmp.endEffectors[i]
                                               .localFrame[j]);
                    V3D v2 =
                        0.2 * ikSolver->kmp.endEffectors[i].targetFrame[j];
                    drawArrow3d(p, v1, 0.01, V3D(1, 0, 0));
                    drawArrow3d(p, v2, 0.01, V3D(0, 1, 0));
                }
            }
        }
    }

    void addJointFKMenu(const char* jName) {
        RobotPose rp(ikSolver->kmp.robot);

        pRBJ j = ikSolver->kmp.robot->getJointByName((char*)jName);
        if (j == NULL)
            return;
        float jointAngle = (float)j->getJointAngle();
        float min = -3.14f;
        float max = 3.14f;
        if (j->hasJointLimits()) {
            min = (float)j->minAngle;
            max = (float)j->maxAngle;
        }

        if (ImGui::SliderFloat(jName, &jointAngle, min, max)) {
            rp.q[ikSolver->kmp.robot->idx(j)] = jointAngle;
            ikSolver->kmp.robot->setPose(rp);
        }
    }

    void addAllJointsToFKMenu(){
        for (auto j : ikSolver->kmp.robot->joints)
            addJointFKMenu(j->name.c_str());
    }

    void addImGuiMenuItems(){
        ImGui::Checkbox("Use Meshes for mouse pick", &checkMeshesForMouseHit);
        ImGui::Checkbox("Use Skeleton for mouse pick", &checkSkeletonForMouseHit);

        if (ImGui::TreeNode("Draw options...")) {
            ikSolver->kmp.robot->df.showDrawOptions();
            ImGui::Separator();
            ImGui::TreePop();
        }

        if (ImGui::TreeNode("Debug options...")) {
            ImGui::Checkbox("Check derivatives", &ikSolver->checkDerivatives);
            ImGui::Checkbox("Print debug info", &ikSolver->minimizer.printOutput);
            ImGui::Checkbox("Test GCRR", &testGCRR);
            ImGui::Separator();
            ImGui::TreePop();
        }

        if (ImGui::TreeNode("IK solver options...")) {
            ImGui::Checkbox("Use 2nd order hessian terms", &ikSolver->kmp.use2ndOrderEETerms);
            ImGui::Checkbox("Freeze root dofs", &ikSolver->kmp.freezeRootDOFs);
            ImGui::Text("Damping regularizer");
            ImGui::SameLine();
            ImGui::InputDouble("#", &ikSolver->minimizer.reg,0, 0, "%.2f");

            ImGui::Text("Obj weight: End Effector Target");
            ImGui::SameLine();
            ImGui::InputDouble("ee_w",&ikSolver->ikObjective.eeObjective.w, 0, 0, "%.2f");
            ImGui::Text("Obj weight: Pose Regularizer");
            ImGui::SameLine();
            ImGui::InputDouble("pose_w",&ikSolver->ikObjective.poseRegularizer.w, 0, 0, "%.2f");
            ImGui::Text("Obj weight: Joint Limits");
            ImGui::SameLine();
            ImGui::InputDouble("jl_w",&ikSolver->ikObjective.jLimitObjective.w, 0, 0, "%.2f");
            ImGui::Separator();
            ImGui::TreePop();
        }

        if (ImGui::TreeNode("Show joint angles")) {
            addAllJointsToFKMenu();
        }
    }

    void drawEETargetWidgets(const Camera& camera){
        if (selectedRB && showTargetSelectionGizmo) {
            ImGui::Begin("Transform Widgets");
            ImGuizmo::BeginFrame();
            V3D scale(1, 1, 1);
            setTransformFromWidgets(
                camera.getViewMatrix(), camera.getProjectionMatrix(), scale,
                targetOrientationForSelectedRB, targetForSelectedPoint);
            if (V3D(oldT, targetForSelectedPoint).norm() > 0.001 ||
                oldQ.angularDistance(targetOrientationForSelectedRB) > 0.001) {
                if (addedEE == false) {
                    Logger::consolePrint("adding end effector\n");
                    addedEE = true;
                    ikSolver->kmp.addEndEffector(EndEffector(selectedRB, selectedRB->getLocalCoordinates(oldT)));
                }
            }

            if (V3D(oldT, targetForSelectedPoint).norm() > 0.001) {
                oldT = targetForSelectedPoint;
                ikSolver->kmp.endEffectors.back().setTargetPosition(
                    targetForSelectedPoint);
                Logger::consolePrint("setting position target\n");
            }
            if (oldQ.angularDistance(targetOrientationForSelectedRB) > 0.001) {
                oldQ = targetOrientationForSelectedRB;
                ikSolver->kmp.endEffectors.back().setTargetOrientation(
                    targetOrientationForSelectedRB);
                Logger::consolePrint("setting orientation target\n");
            }

            ImGui::End();
        }
    }

public:
    IKSolver *ikSolver = nullptr;

    pRB selectedRB = nullptr;
    P3D selectedPoint;
    P3D targetForSelectedPoint, oldT;
    Quaternion targetOrientationForSelectedRB = Quaternion::Identity(),
               oldQ = Quaternion::Identity();
    bool addedEE = false;
    bool showTargetSelectionGizmo = false;
    bool checkMeshesForMouseHit = true;
    bool checkSkeletonForMouseHit = true;
    bool testGCRR = false;

    bool usePersistentEETargets = true;
};
