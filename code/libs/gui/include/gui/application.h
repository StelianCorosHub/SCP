#pragma once

#include <gui/guiElements.h>

#include <gui/camera.h>
#include <gui/inputstate.h>
#include <gui/shader.h>
#include <gui/shadow_casting_light.h>
#include <gui/shadow_map_fbo.h>
#include <utils/logger.h>
#include <utils/timer.h>
#include <imgui.h>
#include <imgui_internal.h>
#include <imgui_widgets/customWidgets.h>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#pragma warning(disable : 4244)

#if defined __APPLE__ && !defined RETINA_SCREEN
#define RETINA_SCREEN
#endif


inline float get_pixel_ratio();

// Sets up a GLFW window, its callbacks and ImGui
class AppInterface {
public:
    GLFWwindow *window;
    int width, height;
    float pixelRatio;

    float clearColor[3] = {0.8f, 0.8f, 0.8f};

    MouseState mouseState;
    KeyboardState keyboardState;

    // FPS params - we will want to display the fps only every 0.3s, otherwise
    // the numbers change way too quickly...
    float avgFPS = 0.0f;
    float percentageOfTimeSpentProcessing = 0.0f;
    float loopRunTime = 0.0f;
    float processRunTime = 0.0f;
    int runningAvgStepCount = 0;

    Timer FPSDisplayTimer;
    Timer tickerPerformanceTimer;
    Timer FPSTimer;

    bool limitFramerate = true;
    int targetFramerate = 60;

    //global app time, advanced with every tick()
    double appTime = 0;
    bool appIsRunning = false;
    bool slowMo = false;

    bool captureScreenshots = false;
    int screenShotCounter = 0;
    char screenshotPath[100] = SCP_DATA_FOLDER "/out/screenshots";

    bool autoManageConsole = true;
    bool showConsole = false;
    int consoleHeight = 250;  // in pixels

    //the two points below will be updated every time the base app draws the app control menu and FPS widget.
    //This will allow other apps to place menus right underneath them, in case that's what they want.
    ImVec2 appMenuBottomLeftCorner;
    ImVec2 FPSLabelBottomLeftCorner;

    void drawFPS();

    void init(const char *title, int width, int height, std::string iconPath);
    void setCallbacks();

    void runMainLoop();

public:
    AppInterface(const char *title, int width, int height, std::string iconPath = SCP_DATA_FOLDER "/crl.png");
    AppInterface(const char *title, std::string iconPath = SCP_DATA_FOLDER "/crl.png");
    virtual ~AppInterface();

    // Called every frame to advance the state of the app forward in time by dt, if the app is running
    virtual void tick(double dt);

    // called once per loop
    virtual void mainLoopStep() {}

    virtual double get_tick_dt() {
        return (1.0 / targetFramerate) * ((slowMo) ? 0.25 : 1.0);
    }

    virtual void draw();
    void swapBuffers();

    void pauseApp(){
        appIsRunning = false;
    }

    void resumeApp(){
        appIsRunning = true;
    }

    virtual void prepareToDraw() {
        glClearColor(clearColor[0], clearColor[1], clearColor[2], 1.00f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT |
                GL_STENCIL_BUFFER_BIT);
    }

    virtual void drawScene() {}

    virtual void drawMenusAndDebugInfo();

    virtual void drawAppSpecificImGuiItems() {}

    //returns the position of the bottomLeftCorner of this window, as we may want to have the app specific menu start right there
    virtual void drawAppControlMenu();

    // override this method to customize placement strategy...
    virtual void drawConsole();
    virtual void resizeWindow(int width, int height);

    virtual void restart() {}

    // return false if the message was not processed, true otherwise
    virtual bool keyPressed(int key, int mods) { return false; }
    virtual bool keyReleased(int key, int mods) { return false; }
    virtual bool mouseButtonPressed(int button, int mods) { return false; }
    virtual bool mouseButtonReleased(int button, int mods) { return false; }
    virtual bool mouseMove(double xpos, double ypos) {return false;}
    virtual bool scrollWheel(double xoffset, double yoffset) {return false;}

    //called when files are dropped onto the app...
    virtual bool drop(int count, const char **filenames) { return false; }

    bool screenshot(const char *filename) const;
};

/**
 * Simple 3d app that has all the basics needed to get started: a camera, a light and a shader
 * */
class BaseApp_Simple : public AppInterface {
public:
    //basic shader that will be used to render the world
    Shader basicShader = Shader(SCP_SHADER_FOLDER "/basic_lighting.vert", SCP_SHADER_FOLDER "/basic_lighting.frag");

    //the shader needs a light...
    Light light = Light();

    //and we'll need a camera too...
    Camera camera;
public:
    BaseApp_Simple(const char *title, int width, int height, std::string iconPath = SCP_DATA_FOLDER "/icons/crl.ico") : AppInterface(title, width, height, iconPath) {
        camera.aspectRatio = float(this->width) / this->height;
    }

    BaseApp_Simple(const char *title, std::string iconPath = SCP_DATA_FOLDER "/icons/crl.ico") : AppInterface(title, iconPath) {
        camera.aspectRatio = float(this->width) / this->height;
    }

    virtual ~BaseApp_Simple() {
    }

    //if nothing upstream handles move and scroll events, pass them to the camera...
    virtual bool mouseMove(double xpos, double ypos) override {
        camera.processMouseMove(mouseState, keyboardState);
        return AppInterface::mouseMove(xpos, ypos);
    }

    virtual bool scrollWheel(double xoffset, double yoffset) override {
        camera.processMouseScroll(xoffset, yoffset);
        return AppInterface::scrollWheel(xoffset, yoffset);
    }

    virtual void resizeWindow(int width, int height) override {
        camera.aspectRatio = float(width) / height;
        return AppInterface::resizeWindow(width, height);
    }

    virtual void prepareToDraw() override{
        AppInterface::prepareToDraw();
        basicShader.setBool("use_textures",false);

        // set up view/projection transformations for the shader - it is up to the app to
        // make sure the camera, light(s) and shader(s) are all compatible with each other...
        basicShader.use();
        basicShader.setMat4("projection", camera.getProjectionMatrix());
        basicShader.setMat4("view", camera.getViewMatrix());
        basicShader.setVec3("camPos", camera.position());
        basicShader.setVec3("lightPos", light.position());
        basicShader.setVec3("lightColor", light.color());

        ShaderFactory::setActiveShader(&basicShader);
    }
};

//perhaps basic crl app with shadows should extend BasicCRLApp?!

class BaseApp : public AppInterface {
public:

    BaseApp(const char *title, int width, int height, std::string iconPath = SCP_DATA_FOLDER "icons/crl.png") : AppInterface(title, width, height, iconPath) {
        if (!shadowMapFBO.Init(this->width, this->height)) {
            std::cout << "Shadow map initialization failed\n";
            exit(0);
        }

        camera.aspectRatio = float(this->width) / this->height;

        camera.distanceToTarget = 8;
        camera.rotAboutUpAxis = 0.75;
        camera.rotAboutRightAxis = 0.5;

        light.s = 0.125f;
        shadowbias = 0.0001f;
    }

    BaseApp(const char *title = "Shadows demo", std::string iconPath = SCP_DATA_FOLDER "/crl.png") : AppInterface(title, iconPath) {

        if (!shadowMapFBO.Init(this->width, this->height)) {
            std::cout << "Shadow map initialization failed\n";
            exit(0);
        }
        camera.aspectRatio = float(this->width) / this->height;

        camera.distanceToTarget = 8;
        camera.rotAboutUpAxis = 0.75;
        camera.rotAboutRightAxis = 0.5;

        light.s = 0.125f;
        shadowbias = 0.0001f;

    }

    //if nothing upstream handles move and scroll events, pass them to the camera...
    virtual bool mouseMove(double xpos, double ypos) override {
        camera.processMouseMove(mouseState, keyboardState);
        return AppInterface::mouseMove(xpos, ypos);
    }

    virtual bool scrollWheel(double xoffset, double yoffset) override {
        camera.processMouseScroll(xoffset, yoffset);
        return AppInterface::scrollWheel(xoffset, yoffset);
    }

    virtual void resizeWindow(int width, int height) override {
        if (!shadowMapFBO.Init(this->width, this->height)) {
            std::cout << "Shadow map initialization failed\n";
            exit(0);
        }

        camera.aspectRatio = float(width) / height;
        return AppInterface::resizeWindow(width, height);
    }

    virtual void drawScene() override {
        shadowPass();
        renderPass();
    }

    virtual void prepareToDraw() override{
        AppInterface::prepareToDraw();
    }

    void fixLightsAndCam(const P3D& p, bool adjustCameraTarget){
        //make the camera and light focus in on the robot...
        //light does by default, to get the shadows to behave...
        light.target.x() = p.x;
        light.target.z() = p.z;
        if (adjustCameraTarget){
            camera.target.x = p.x;
            camera.target.z = p.z;
        }
    }

    void shadowPass() {
        shadowMapFBO.BindForWriting();
        glClear(GL_DEPTH_BUFFER_BIT);
        shadowMapRenderer.use();
        shadowMapRenderer.setMat4("projection",light.getOrthoProjectionMatrix());
        shadowMapRenderer.setMat4("view", light.getViewMatrix());
        glViewport(0, 0, shadowMapFBO.bufferWidth, shadowMapFBO.bufferHeight);

        ShaderFactory::setActiveShader(&shadowMapRenderer);
        drawShadowCastingObjects();

        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glUseProgram(0);
#ifdef RETINA_SCREEN
        // temporal one... need more investigation
        glViewport(0, 0, width * 2, height * 2);
#else
        glViewport(0, 0, width, height);
#endif
    }

    void renderPass() {
        shadowMapFBO.BindForReading(GL_TEXTURE0);

#define SETUP_SHADER(shader)                                        \
    shader.use();                                                   \
    shader.setMat4("projection", camera.getProjectionMatrix());     \
    shader.setMat4("view", camera.getViewMatrix());                 \
    shader.setVec3("camPos", camera.position());                    \
    (shader).setVec3("lightPos", light.position());                 \
    shader.setVec3("lightColor", light.color());

        // set up shaders
        SETUP_SHADER(shadowShader);
        shadowShader.setMat4("lightProjection", light.getOrthoProjectionMatrix());
        shadowShader.setMat4("lightView", light.getViewMatrix());
        shadowShader.setInt("shadowMap", 0);
        shadowShader.setFloat("bias", shadowbias);

        SETUP_SHADER(basicShader);
        // better lighting approximation here so that regions of the model do
        // not remain forever shaded dark...
        basicShader.setVec3("lightPos", camera.position());

        ShaderFactory::setActiveShader(&shadowShader);
        drawObjectsWithShadows();
        ShaderFactory::setActiveShader(&basicShader);
        drawObjectsWithoutShadows();

    }

    // objects drawn in this function will cast a shadow
    virtual void drawShadowCastingObjects() = 0;

    // objects drawn here will have shadows cast on them
    virtual void drawObjectsWithShadows() = 0;

    // objects drawn here will not have shadows cast on them
    virtual void drawObjectsWithoutShadows() = 0;

public:
    //basic shader that will be used to render the world
    Shader basicShader = Shader(SCP_SHADER_FOLDER "/basic_lighting.vert", SCP_SHADER_FOLDER "/basic_lighting.frag");

    //the shader needs a light...
    ShadowCastingLight light = ShadowCastingLight();

    //and we'll need a camera too...
    Camera camera;

    ShadowMapFBO shadowMapFBO;
    float shadowbias = 0.0001f;

    Shader shadowShader =
        Shader(SCP_SHADER_FOLDER "/basic_lighting.vert",
               SCP_SHADER_FOLDER "/basic_shadow_lighting.frag");
    Shader shadowMapRenderer = Shader(SCP_SHADER_FOLDER "/basic_lighting.vert",
                                      SCP_SHADER_FOLDER "/render_shadow.frag");

};

