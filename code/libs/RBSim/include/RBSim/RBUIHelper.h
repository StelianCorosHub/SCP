#pragma once


/*
 * TODO:
 * - pick body parts, tell the name, the parent, the child joints, the mass...
 * - select joints as well - tell the parent and the child
 * - pick the meshes and the joints. Do we go back to picking skeleton?!
*/


/**
 * UI helper class for interactions with rigid body systems.
 *
 * Left Click + drag: set and move target on selected robot body part
 */
class RBUIHelper {
public:
    RBUIHelper(RBSimEngine* rbEngine) {
        this->simEngine = rbEngine;
    }

    RBUIHelper(){
        this->simEngine = nullptr;
    }

    ~RBUIHelper(){}

    bool mouseDrag(MouseState mouseState, const Ray& mouseRay) {
        if (simEngine == nullptr) return false;

        if (mouseState.lButtonPressed == true) {
            if (selectedRB != NULL) {
                mouseRay.getDistanceToPoint(selectedRB->getWorldCoordinates(selectedPoint), &targetForSelectedPoint);
                return selectedRB != NULL;
            }
        }
        return false;
    }

    bool mouseButtonReleased(int button, int mods) {
        if (simEngine == nullptr) return false;

        if (button == GLFW_MOUSE_BUTTON_LEFT && selectedRB) {
            selectedRB->selected = false;
            selectedRB = nullptr;
        }
        return true;
    }

    bool mouseButtonPressed(int button, int mods, const Ray& mouseRay) {
        if (simEngine == nullptr) return false;

        if (button == GLFW_MOUSE_BUTTON_LEFT || button == GLFW_MOUSE_BUTTON_RIGHT) {
            if (selectedRB != NULL)
                selectedRB->selected = false;

            selectedRB = simEngine->getFirstRBHitByRay(mouseRay, selectedPoint);

            if (selectedRB != NULL) {

                selectedRB->selected = true;
                targetForSelectedPoint = selectedPoint;
                selectedPoint = selectedRB->getLocalCoordinates(selectedPoint);

                Logger::consolePrint("clicked on %s at %lf %lf %lf\n", selectedRB->name.c_str(),
                                     selectedPoint.x,
                                     selectedPoint.y,
                                     selectedPoint.z);
            }

            return true;
        }
        return false;
    }

    // objects drawn here will not have shadows cast on them
    void drawDebugInfo() {
        if (selectedRB) {
            drawSphere(targetForSelectedPoint, 0.05f);
            drawCapsule(selectedRB->getWorldCoordinates(selectedPoint), targetForSelectedPoint, 0.01f);
        }
    }

    void addImGuiMenuItems(){
        ImGui::Checkbox("Use Meshes for mouse pick", &checkMeshesForMouseHit);
        ImGui::Checkbox("Use Skeleton for mouse pick", &checkSkeletonForMouseHit);

        ImGui::InputDouble("InteractionForce kp:", &kp);
        ImGui::InputDouble("InteractionForce kd:", &kd);
    }

    V3D getInteractionForce() {
        if (selectedRB == nullptr)
            return V3D();

        return -kp * V3D(targetForSelectedPoint, selectedRB->getWorldCoordinates(selectedPoint)) - kd * selectedRB->getVelocityForLocalCoordsPoint(selectedPoint);
    }

    void applyInteractionForce() {
        if (simEngine == nullptr || selectedRB == nullptr)
            return;
        simEngine->addForceToRB(selectedRB, getInteractionForce(), selectedPoint);
    }

public:
    RBSimEngine *simEngine = nullptr;
    pRB selectedRB = nullptr;

    P3D selectedPoint;
    P3D targetForSelectedPoint;

    bool checkMeshesForMouseHit = true;
    bool checkSkeletonForMouseHit = true;

    double kp = 1000;
    double kd = 10;
};

