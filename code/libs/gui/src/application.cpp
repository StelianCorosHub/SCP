#include <gui/application.h>
#include <imgui_widgets/implot.h>

// defined in model.cpp
//#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#include <gui/model.h>

float get_pixel_ratio() {
#ifdef RETINA_SCREEN
    return 1.f;
#endif
    GLFWmonitor *monitor = glfwGetPrimaryMonitor();
    if (monitor == nullptr) throw "Primary monitor not found.";
    float xscale, yscale;
    glfwGetMonitorContentScale(monitor, &xscale, &yscale);
    return xscale;
}

void AppInterface::init(const char *title, int width, int height,
                    std::string iconPath) {
    // glfw: initialize and configure

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SAMPLES, 8);

    glfwWindowHint(GLFW_FOCUSED, GLFW_TRUE);

#ifdef SINGLE_BUFFER
    glfwWindowHint(GLFW_DOUBLEBUFFER, GL_FALSE);  // turn off framerate limit
#endif

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT,
                   GL_TRUE);  // fix compilation on OS X
#endif

    // get pixel ratio and adjust width/height
    pixelRatio = get_pixel_ratio();

    // glfw window creation
    window = glfwCreateWindow(width, height, title, nullptr, nullptr);
    if (window == nullptr) {
        glfwTerminate();
        throw std::runtime_error("Failed to create GLFW window");
    }
    glfwMakeContextCurrent(window);

    glfwGetWindowSize(window, &this->width, &this->height);

    // app icon
    if (iconPath != "") {
        GLFWimage image;
        image.pixels = stbi_load(iconPath.c_str(), &image.width, &image.height,
                                 nullptr, 4);
        glfwSetWindowIcon(window, 1, &image);
        stbi_image_free(image.pixels);
    }

    // glad: load all OpenGL function pointers
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        throw std::runtime_error("Failed to initialize GLAD");
    }

    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
    glEnable(GL_MULTISAMPLE);

    // Setup Dear ImGui binding
    const char *glsl_version = "#version 150";
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGuiIO &io = ImGui::GetIO();
    ImFontConfig cfg;
    cfg.SizePixels = 40 * pixelRatio;

    ImFont *imFont = io.Fonts->AddFontFromFileTTF(
        SCP_IMGUI_FONT_FOLDER "/Roboto-Medium.ttf", 15.0f * pixelRatio, &cfg);
    imFont->DisplayOffset.y = pixelRatio;
    ImGuiStyle &style = ImGui::GetStyle();
    style.ScaleAllSizes(pixelRatio);

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);

    setCallbacks();

    glEnable(GL_DEPTH_TEST);

    //create a dummy first texture to remove shader warnings
    unsigned int textureID;
    glGenTextures(1, &textureID);
    glBindTexture(GL_TEXTURE_2D, textureID);

    unsigned char ch = 0;
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, 1, 1, 0, GL_RED, GL_UNSIGNED_BYTE, &ch);

    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
}

AppInterface::AppInterface(const char *title, int width, int height,
                   std::string iconPath) {
    if (!glfwInit()) {
        // An error occured
        std::cout << "GLFW initialization failed\n";
        exit(0);
    }

    init(title, width, height, iconPath);
}

AppInterface::AppInterface(const char *title, std::string iconPath) {

    if (!glfwInit()) {
        // An error occured
        std::cout << "GLFW initialization failed\n";
        exit(0);
    }

    const GLFWvidmode *mode = glfwGetVideoMode(glfwGetPrimaryMonitor());

#ifdef __APPLE__
    int borderLeft = 0;
    int borderTop = 0;
    int borderRight = 0;
    int borderBottom = 0;
#else
    int borderLeft = 2;
    int borderTop = 70;
    int borderRight = 2;
    int borderBottom = 105;
#endif

    init(title, (mode->width - borderLeft - borderRight),
         (mode->height - borderTop - borderBottom), iconPath);
    glfwSetWindowPos(window, borderLeft, borderTop);

}

AppInterface::~AppInterface() {
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
}

void AppInterface::setCallbacks() {
    glfwSetWindowUserPointer(window, this);

    glfwSetErrorCallback([](int error, const char *description) {
        std::cout << "Error " << error << ": " << description << std::endl;
    });

    glfwSetFramebufferSizeCallback(window, [](GLFWwindow *window, int width,
                                              int height) {
        auto app = static_cast<AppInterface *>(glfwGetWindowUserPointer(window));
        app->resizeWindow(width, height);
        //by default app is kinda blocked while resizing, which looks bad... so force a redraw here.
        app->draw();
        app->swapBuffers();
    });

    glfwSetKeyCallback(window, [](GLFWwindow *window, int key, int scancode,
                                  int action, int mods) {
        auto app = static_cast<AppInterface *>(glfwGetWindowUserPointer(window));
        app->keyboardState[key] = (action != GLFW_RELEASE);

        ImGui_ImplGlfw_KeyCallback(window, key, scancode, action, mods);
        if (ImGui::GetIO().WantCaptureKeyboard || ImGui::GetIO().WantTextInput) {
            return;
        }

        if (key == GLFW_KEY_ESCAPE) {
            glfwSetWindowShouldClose(window, GL_TRUE);
            return;
        }

        if (key == GLFW_KEY_SPACE && action == GLFW_PRESS) {
            app->appIsRunning = !app->appIsRunning;
            return;
        }

        if (key == GLFW_KEY_GRAVE_ACCENT && action == GLFW_PRESS) {
            app->showConsole = !app->showConsole;
            return;
        }

        if (action == GLFW_PRESS) app->keyPressed(key, mods);

        if (action == GLFW_RELEASE) app->keyReleased(key, mods);
    });

    glfwSetMouseButtonCallback(window, [](GLFWwindow *window, int button,
                                          int action, int mods) {
        if (action == GLFW_PRESS || action == GLFW_RELEASE) {
            glfwFocusWindow(window);
        }

        double xPos, yPos;
        glfwGetCursorPos(window, &xPos, &yPos);
#ifdef RETINA_SCREEN
        xPos *= 2;
        yPos *= 2;
#endif
        auto app = static_cast<AppInterface *>(glfwGetWindowUserPointer(window));
        app->mouseState.onMouseClick(xPos, yPos, button, action, mods);

        if (ImGui::GetIO().WantCaptureMouse) {
            ImGui_ImplGlfw_MouseButtonCallback(window, button, action, mods);
            return;
        }

        if (action == GLFW_PRESS) app->mouseButtonPressed(button, mods);

        if (action == GLFW_RELEASE) app->mouseButtonReleased(button, mods);
    });

    glfwSetCursorPosCallback(window, [](GLFWwindow *window, double xpos,
                                        double ypos) {
        auto app = static_cast<AppInterface *>(glfwGetWindowUserPointer(window));
#ifdef RETINA_SCREEN
        xpos *= 2;
        ypos *= 2;
#endif
        app->mouseState.onMouseMove(xpos, ypos);

        if (ImGui::GetIO().WantCaptureMouse) return;

        app->mouseMove(xpos, ypos);
    });

    glfwSetScrollCallback(window, [](GLFWwindow *window, double xoffset,
                                     double yoffset) {
        if (ImGui::GetIO().WantCaptureMouse) {
            ImGui_ImplGlfw_ScrollCallback(window, xoffset, yoffset);
            return;
        }

        auto app = static_cast<AppInterface *>(glfwGetWindowUserPointer(window));
        app->scrollWheel(xoffset, yoffset);
    });

    glfwSetDropCallback(window, [](GLFWwindow *window, int count,
                                   const char **filenames) {
        auto app = static_cast<AppInterface *>(glfwGetWindowUserPointer(window));
        app->drop(count, filenames);
    });
}

void AppInterface::runMainLoop() {
    glfwSwapInterval(0);  //disable waiting for framerate of glfw window
    //make this window active, for some reason at times it likes to lose focus before loading is complete, at least on macos
    glfwFocusWindow(window);

    while (!glfwWindowShouldClose(window)) {
        if (FPSDisplayTimer.timeEllapsed() > 0.33) {
            FPSDisplayTimer.restart();
            if (runningAvgStepCount > 0) {
                avgFPS = 1.0 / (loopRunTime / runningAvgStepCount);
                percentageOfTimeSpentProcessing = processRunTime / loopRunTime;
            } else {
                avgFPS = -1;
                percentageOfTimeSpentProcessing = -1;
            }
            loopRunTime = 0;
            processRunTime = 0;
            runningAvgStepCount = 0;
        }

        runningAvgStepCount++;
        loopRunTime += FPSTimer.timeEllapsed();
        FPSTimer.restart();

        //TODO: we might want a way here to make sure the entire app runs ok even with variable frame rates,
        // a mechanism that aims to keep appTime synchronized with wall-clock time.

        mainLoopStep();

        if (appIsRunning) {
            tickerPerformanceTimer.restart();
            tick(get_tick_dt());
            appTime += get_tick_dt();
            processRunTime += tickerPerformanceTimer.timeEllapsed();
        }

        // draw, swap buffers and poll IO events (keys pressed/released, mouse moved etc.)

        draw();
        swapBuffers();
        glfwPollEvents();

        if (limitFramerate)
            while (FPSTimer.timeEllapsed() < (1.0f / targetFramerate));

        if (captureScreenshots) {
            char filename[1000];
            sprintf(filename, "%s_%04d.png", screenshotPath, screenShotCounter);
            screenshot(filename);
            screenShotCounter++;
        }
    }

    // glfw: terminate, clearing all previously allocated GLFW resources.
    glfwTerminate();
}

// Called every frame to advance the state of the app forward in time by dt
void AppInterface::tick(double dt) {}

void AppInterface::drawAppControlMenu() {
    ImGui::SetNextWindowPos(ImVec2(0, 0), ImGuiCond_Once);
    ImGui::Begin("App Control Menu");

    ImGui::Text("Play Menu:");
    ImGui::SameLine();

    if (DoubleBackButton("Reload")) {
        restart();
    }
    ImGui::SameLine();
    PlayPauseButton("PlayButton", &appIsRunning);

    if (appIsRunning == false) {
        ImGui::SameLine();
        if (ImGui::ArrowButton("tmpFF", ImGuiDir_Right)) {
            tick(get_tick_dt());
            appTime += get_tick_dt();
        }
    }

    ImGui::Checkbox("SlowMo Mode", &slowMo);
    ImGui::Checkbox("Capture Screenshots", &captureScreenshots);

    ImGui::Checkbox("Draw Console", &showConsole);
    if (showConsole) ImGui::Checkbox("Auto-manage Console", &autoManageConsole);

    appMenuBottomLeftCorner.x = ImGui::GetWindowPos().x;
    appMenuBottomLeftCorner.y =
        ImGui::GetWindowPos().y + ImGui::GetWindowHeight();

    ImGui::End();
}

void AppInterface::draw() {
    prepareToDraw();

    drawScene();

    drawMenusAndDebugInfo();
}

void AppInterface::swapBuffers() {
#ifdef SINGLE_BUFFER
    glFlush();
#else
    glfwSwapBuffers(window);
#endif
}

void AppInterface::drawMenusAndDebugInfo() {
    // draw ImGui menu
    {
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        drawFPS();

        drawAppControlMenu();

        drawConsole();

        drawAppSpecificImGuiItems();

        ImGui::EndFrame();
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    }
}

void AppInterface::drawFPS() {
    ImGui::SetNextWindowPos(ImVec2(this->width - pixelRatio * 320, 0),
                            ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(pixelRatio * 320, pixelRatio * 80),
                             ImGuiCond_Always);
    ImGui::SetNextWindowCollapsed(true, ImGuiCond_Once);
    char title[100];
    sprintf(title, "FPS: %.2f###FPS", avgFPS);
    ImGui::Begin(title);

    if (appIsRunning)
        ImGui::Text("Time spent processing: %.2f%%",
                    percentageOfTimeSpentProcessing);

    ImGui::Checkbox("Limit FPS", &limitFramerate);
    ImGui::SameLine(pixelRatio * 100);
    if (limitFramerate) ImGui::InputInt("", &targetFramerate);

    FPSLabelBottomLeftCorner.x = ImGui::GetWindowPos().x;
    FPSLabelBottomLeftCorner.y =
        ImGui::GetWindowPos().y + ImGui::GetWindowHeight();

    ImGui::End();
}

void AppInterface::drawConsole() {
    if (showConsole == false) return;

    if (autoManageConsole == false) {
        ImGui::SetNextWindowSize(ImVec2(this->width, pixelRatio * 335),
                                 ImGuiCond_Once);
        ImGui::SetNextWindowPos(ImVec2(0, this->height - pixelRatio * 350),
                                ImGuiCond_Once);
    }
    ImGui::Begin("Console");
    if (autoManageConsole == true) {
        if (ImGui::IsWindowCollapsed()) {
            ImGui::SetWindowPos(ImVec2(this->width - pixelRatio * 300,
                                       this->height - pixelRatio * 20),
                                ImGuiCond_Always);
            ImGui::SetWindowSize(ImVec2(pixelRatio * 300, pixelRatio * 80),
                                 ImGuiCond_Always);
        } else {
            ImGui::SetWindowPos(
                ImVec2(0, this->height - pixelRatio * consoleHeight),
                ImGuiCond_Always);
            ImGui::SetWindowSize(
                ImVec2(this->width, pixelRatio * consoleHeight),
                ImGuiCond_Always);
        }
    }

    for (const Logger::ConsoleText &cText : Logger::consoleOutput) {
        ImVec4 color(cText.color.x(), cText.color.y(), cText.color.z(), 1.0);
        ImGui::TextColored(color, "%s", cText.text.c_str());
    }
    ImGui::End();
}

void AppInterface::resizeWindow(int width, int height) {
    this->width = width;
    this->height = height;

#ifdef RETINA_SCREEN
    this->width /= 2.0;
    this->height /= 2.0;
#endif

    glViewport(0, 0, width, height);
}

bool AppInterface::screenshot(const char *filename) const {
    std::vector<unsigned char> pixels(width * height * 3);
    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, &pixels[0]);
    stbi_flip_vertically_on_write(1);
    return (bool)stbi_write_png(filename, width, height, 3, &pixels[0], 0);
}
