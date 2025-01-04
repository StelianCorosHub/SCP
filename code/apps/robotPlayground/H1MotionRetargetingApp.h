#pragma once

#include <gui/application.h>

#include <utils/BVHReader.h>

#include "H1.h"
#include "H1MotionAnalyzer.h"

/**
 * App that retargets BVH files (from the 100styles dataset) to full-body motions for H1.
 **/
class H1MotionRetargetingApp : public BaseApp {
public:
    H1MotionRetargetingApp(const char* title = "Motion Retargetting App") : BaseApp(title) {

//        loadMocapClip(SCP_DATA_FOLDER "/bvh/human/100styles/neutral/Transitions.bvh");
        loadMocapClip(SCP_DATA_FOLDER "/bvh/human/100styles/neutral/Forwards Walking.bvh");
//        loadMocapClip(SCP_DATA_FOLDER "/bvh/human/100styles/neutral/Backwards Walking.bvh");
//        loadMocapClip(SCP_DATA_FOLDER "/bvh/human/100styles/neutral/Sidestep Walking.bvh");
//        loadMocapClip(SCP_DATA_FOLDER "/bvh/human/100styles/neutral/Forwards Running.bvh");

//        loadMocapClip(SCP_DATA_FOLDER "/bvh/human/100styles/allStyles/Tiptoe/Tiptoe_FW.bvh");

        START_FIDX = 0;
        END_FIDX = mocapClip.getFrameCount();

        fp = fopen(SCP_DATA_FOLDER "/out/H1Walk", "w");
        fprintf(fp, "# Root position x, y, z\n# Root orientation w, x, y, z\n");
        for (auto j : h1->joints)
            fprintf(fp, "# %s\n", j->name.c_str());

        updateMocapPose();
    }

    void loadMocapClip(const char* fName) {
        Logger::consolePrint("Bip has a mass of %lfkg and a height of %2.2lfm\n", h1->getTotalMass(), h1->getHeight());
        mocapClip.loadFromFile(fName);
        Logger::consolePrint("Loaded %d frames for a skeleton that is %2.2lfm tall and has %d degrees of freedom. Motion data is available at %2.2lfHz.\n", mocapClip.motionData.size(), mocapClip.skeleton.getUnscaledHeight(), mocapClip.skeleton.getChannelCount(), 1.0 / mocapClip.frame_dt);

        ma.initialize(h1, &mocapClip);
    }

    virtual void updateMocapPose() {
        clamp(fIdx, START_FIDX, END_FIDX - 1);
        msp = (showZeroPose) ? ma.getDefaultMocapPose() : ma.getMocapPose(fIdx);
        ma.setH1PoseFromMocapSkeletonPose(msp);

        fixLightsAndCam(msp.getRootPosition(), trackSkeleton);
    }

    virtual bool keyReleased(int key, int mods) override {
        if (key == GLFW_KEY_LEFT){
            fIdx -= 2;
            tick(get_tick_dt());
            return true;
        }
        if (key == GLFW_KEY_RIGHT){
            tick(get_tick_dt());
            return true;
        }
        return BaseApp::keyReleased(key, mods);
    }

    virtual void restart() override {
        fIdx = START_FIDX;
        updateMocapPose();
    }

    // Called every frame to advance the state of the app forward in time by dt
    void tick(double dt) override {
        fIdx++;
        updateMocapPose();

        if (recordRobotPoses) {
            fprintf(fp, "\n\n");
            fprintf(fp, "%lf %lf %lf\n%lf %lf %lf %lf\n", h1->root->position.x, h1->root->position.y, h1->root->position.z,
                    h1->root->orientation.w(), h1->root->orientation.x(), h1->root->orientation.y(), h1->root->orientation.z());
            for (auto j : h1->joints)
                fprintf(fp, "%lf\n", j->getJointAngle());
            fflush(fp);
        }
    }

    virtual ~H1MotionRetargetingApp() override {}

    // objects drawn here will not have shadows cast on them
    virtual void drawObjectsWithoutShadows() override {
        if (showH1)
            h1->draw();
        if (showMocapSkeleton)
            drawBVHPose(msp, 0.025 * 1.8, 0.4, selectedBoneIndex);

    }

    // objects drawn in this function will cast a shadow
    virtual void drawShadowCastingObjects() override {
        if (showMocapSkeleton)
            drawBVHPose(msp, 0.025 * 1.8);
        if (showH1)
            h1->draw();
    }

    // objects drawn here will have shadows cast on them
    virtual void drawObjectsWithShadows() override {
        ground.draw();
    }

    virtual void drawAppSpecificImGuiItems() override {
        ImGui::SetNextWindowPos( appMenuBottomLeftCorner, ImGuiCond_Once);

        ImGui::Begin("Main App Menu");

        if (ToggleButton("camtrack", &trackSkeleton))
            if (trackSkeleton)
                fixLightsAndCam(msp.getRootPosition(), trackSkeleton);
        ImGui::SameLine();
        ImGui::Text("Track skeleton");

        if (ImGui::TreeNode("Draw options...")) {
            h1->df.showDrawOptions();
            ImGui::TreePop();
            ImGui::Separator();
        }

        if (ToggleButton("defaultpose", &showZeroPose))
            updateMocapPose();
        ImGui::SameLine(); ImGui::Text("Show Zero Pose");

        ToggleButton("showh1", &showH1); ImGui::SameLine();
        ImGui::Text("Show Robot");

        ToggleButton("showmmocap", &showMocapSkeleton); ImGui::SameLine();
        ImGui::Text("Show Mocap Skeleton");

        ToggleButton("recordposes", &recordRobotPoses); ImGui::SameLine();
        ImGui::Text("RecordMocapPoses");

        if (ImGui::Button("<<")){
            fIdx -= 6;
            tick(get_tick_dt());
        } ImGui::SameLine();
        if (ImGui::Button("<")){
            fIdx -= 2;
            tick(get_tick_dt());
        } ImGui::SameLine();
        if (ImGui::Button(">")){
            fIdx += 0;
            tick(get_tick_dt());
        } ImGui::SameLine();
        if (ImGui::Button(">>")){
            fIdx += 4;
            tick(get_tick_dt());
        }

        if (ImGui::SliderInt("Frame Idx", &fIdx, START_FIDX, END_FIDX - 1))
            updateMocapPose();

        ImGui::SliderInt("Selected Bone Idx", &selectedBoneIndex, -1, msp.getBoneCount()-1);

        ImGui::End();
    }

public:
    SimpleGroundModel ground;
    pH1 h1 = pH1(new H1());
    H1MotionAnalyzer ma;

    BVHClip mocapClip;

    int fIdx = 0;

    bool showZeroPose = false;

    bool trackSkeleton = true;

    bool showMocapSkeleton = true;
    bool showH1 = true;

    bool recordRobotPoses = false;
    FILE* fp = nullptr;

    MocapSkeletonPose msp;

    int selectedBoneIndex = -1;

    int START_FIDX = 0, END_FIDX = 0;
};
