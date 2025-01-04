#pragma once

#include <gui/application.h>
#include <RBSim/ArticulatedFigure.h>
#include <control/GeneralizedCoordinatesRepresentation.h>

/**
 * Basic robot viewer
 **/
class RobotViewer : public BaseApp {
public:
    RobotViewer(const char* title = "Robot Viewer") : BaseApp(title) {
        camera.distanceToTarget = 8;
        camera.rotAboutUpAxis = 0.75;
        camera.rotAboutRightAxis = 0.5;

        light.s = 0.125f;
        shadowbias = 0.0001f;

        showConsole = true;

        loadRobot(SCP_DATA_FOLDER "/Robots/h1/h1.rbs");
    }

    void loadRobot(const char* fName){
       robot = pRobot(new Robot(fName));
       robot->setDefaultState();
       if (robotName != std::string() + fName && robotName.c_str() != fName)
           robotName = std::string() + fName;
       Logger::consolePrint("Robot has a mass of %lfkg, and it consists of %d body parts and %d joints\n", robot->getTotalMass(), robot->getRBCount(), robot->getJointCount());
    }

    virtual ~RobotViewer() override {}

    virtual void restart() override {
        loadRobot(robotName.c_str());
    }

    bool mouseButtonPressed(int button, int mods) override {
        if (selectedRB.get() != nullptr)
            selectedRB->selected = false;

        P3D p;
        selectedRB = robot->getFirstRBHitByRay(camera.getRayFromScreenCoordinates(mouseState.lastMouseX, mouseState.lastMouseY), p);

        if (selectedRB.get() != nullptr) {
            selectedRB->selected = true;
            p = selectedRB->getLocalCoordinates(p);
            Logger::consolePrint("Clicked on '%s' - COM position: %lf %lf %lf, local coords: %lf %lf %lf\n", selectedRB->name.c_str(),
                                 selectedRB->position.x, selectedRB->position.y, selectedRB->position.z, p.x, p.y, p.z);
        }

        return BaseApp::mouseButtonPressed(button, mods);
    }

    virtual void drawAppSpecificImGuiItems() override {
        ImGui::SetNextWindowPos( appMenuBottomLeftCorner, ImGuiCond_Once);

        ImGui::Begin("Main App Menu");

        if (ImGui::Button("Test GeneralizedCoordinatedRepresentation"))
            testGeneralizedCoordinateRepresentation(robot.get());

        if (ImGui::TreeNode("Draw options...")) {
            robot->df.showDrawOptions();
            ImGui::TreePop();
            ImGui::Separator();
        }

        if (ImGui::Button("Set Zero Pose"))
            robot->setZeroState();
        if (ImGui::Button("Set DefaultPose"))
            robot->setDefaultState();
        if (ImGui::Button("Set RandomPose")) {
            robot->setRandomState(0.1);
        }

        if (ImGui::TreeNode("Joint Angles")) {
            bool dirty = false;
            RobotState rs = robot->getState();
            for (int i = 0; i < robot->getJointCount(); i++){
                float jointAngle = (float) rs.q[i];
                float min = (float) -PI;
                float max = (float) PI;
                if (robot->joints[i]->hasJointLimits()) {
                    min = (float) robot->joints[i]->minAngle;
                    max = (float) robot->joints[i]->maxAngle;
                }

                if (ImGui::SliderFloat(robot->joints[i]->name.c_str(), &jointAngle, min, max)){
                    dirty = true;
                    rs.q[i] = jointAngle;
                }
            }

            if (dirty)
                robot->setState(rs);

            ImGui::TreePop();
            ImGui::Separator();
        }

        ImGui::End();
    }

    // objects drawn in this function will cast a shadow
    virtual void drawShadowCastingObjects() override {
        robot->draw();
    }

    // objects drawn here will have shadows cast on them
    virtual void drawObjectsWithShadows() override {
        ground.draw();
    }

    // objects drawn here will not have shadows cast on them
    virtual void drawObjectsWithoutShadows() override {
        robot->draw();
    }

    virtual bool drop(int count, const char** fileNames) override {
        loadRobot(fileNames[count - 1]);
        return true;
    }

public:
    SimpleGroundModel ground;

    pRobot robot;
    std::string robotName;

    pRB selectedRB;

};


