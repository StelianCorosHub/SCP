#include <gui/model.h>

#include <iostream>

#include "RobotViewer.h"
#include "H1IKApp.h"
#include "H1MotionRetargetingApp.h"

int main() {
//    RobotViewer app;
//    H1IKApp app;

    H1MotionRetargetingApp app;


    app.runMainLoop();

    return 0;
}

