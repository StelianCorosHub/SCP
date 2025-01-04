#pragma once

#include <gui/application.h>
#include <gui/camera.h>
#include <gui/light.h>
#include <gui/renderer.h>
#include <gui/shader.h>
#include <control/kinematics/IKSolver.h>
#include <RBSim/ArticulatedFigure.h>
#include <utils/logger.h>

#include <control/kinematics/IKUIHelper.h>

#include "H1.h"

/**
 * H1 IK demo app.
 **/
class H1IKApp : public BaseApp {
public:
    H1IKApp(const char *title = "SCP - H1 IK App") : BaseApp(title) {
        camera.distanceToTarget = 5;
        camera.rotAboutUpAxis = 0.75;
        camera.rotAboutRightAxis = 0.5;

        light.s = 0.03f;
        shadowbias = 0.0001f;

        h1->setDefaultState();

        ikSolver = new IKSolver(h1.get(), true);
        uiHelper = IKUIHelper(ikSolver);
    }

    virtual ~H1IKApp() override { delete ikSolver; }

    virtual void restart() override{
        h1->setDefaultState();
    }

    bool mouseMove(double xpos, double ypos) override {
        if (uiHelper.mouseDrag(mouseState, camera.getRayFromScreenCoordinates(xpos, ypos)))
            return true;

        return BaseApp::mouseMove(xpos, ypos);
    }

    bool mouseButtonReleased(int button, int mods) override {
        if (uiHelper.mouseButtonReleased(button, mods))
            return true;
        return BaseApp::mouseButtonReleased(button, mods);
    }

    bool mouseButtonPressed(int button, int mods) override {
        if (uiHelper.mouseButtonPressed(button, mods, camera.getRayFromScreenCoordinates(mouseState.lastMouseX, mouseState.lastMouseY)))
            return true;
        return BaseApp::mouseButtonPressed(button, mods);
    }

    void tick(double dt) override {
        ikSolver->solve();
        if (uiHelper.testGCRR)
            testGeneralizedCoordinateRepresentation(h1.get());
    }

    // objects drawn in this function will cast a shadow
    virtual void drawShadowCastingObjects() override {
        h1->draw();
    }

    // objects drawn here will have shadows cast on them
    virtual void drawObjectsWithShadows() override {
//        robot.draw();
        if (drawGround)
            ground.draw(V3D(0.9, 0.9, 0.9));
    }

    // objects drawn here will not have shadows cast on them
    virtual void drawObjectsWithoutShadows() override {
        h1->draw();
        uiHelper.drawDebugInfo();
    }

    virtual void drawAppSpecificImGuiItems() override {
        ImGui::SetNextWindowPos( appMenuBottomLeftCorner, ImGuiCond_Once);

        ImGui::Begin("Main Menu");
        ImGui::Checkbox("Draw Ground", &drawGround);
        uiHelper.addImGuiMenuItems();

        ImGui::End();

        uiHelper.drawEETargetWidgets(camera);
    }

public:
    SimpleGroundModel ground;
    bool drawGround = true;

    IKSolver *ikSolver = nullptr;
    IKUIHelper uiHelper;
    pH1 h1 = pH1(new H1());
};
