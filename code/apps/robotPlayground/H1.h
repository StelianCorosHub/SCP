#pragma once

#include <RBSim/ArticulatedFigure.h>
#include <RBSim/RBSLoader.h>
#include <utils/spline.h>

#define H1_FORWARD V3D(0,0,1)
#define H1_RIGHT V3D(-1,0,0)
#define H1_UP UP

/**
 * This class implements specific methods for H1, a humanoid robot.
 */
class H1 : public ArticulatedFigure {

public:
    /** the constructor */
    H1(bool loadVisuals = true){
        RBSLoader* rbsLoader = new RBSLoader(SCP_DATA_FOLDER "/robots/h1/h1.rbs", loadVisuals);
        rbsLoader->populate(this);
        delete rbsLoader;
    }

    /** the destructor */
    virtual ~H1(void){}

    double getHeight(){
        RobotState s(this);

        setZeroState();

        pRB head = getRBByName("torso_link");
        pRB rFoot = getRBByName("right_ankle_roll_link");

        P3D headM = head->getWorldCoordinates(head->endEffectorPoints[0].pos);
        P3D toesM = rFoot->getWorldCoordinates(rFoot->endEffectorPoints[0].pos);

        setState(s);
        return V3D(toesM, headM).dot(UP);
    }

    P3D getRightHeelPosition() {
        pRB foot = getRBByName("right_ankle_roll_link");
        return (foot->getWorldCoordinates(foot->endEffectorPoints[0].pos) + foot->getWorldCoordinates(foot->endEffectorPoints[1].pos)) / 2.0;
    }

    P3D getLeftHeelPosition() {
        pRB foot = getRBByName("left_ankle_roll_link");
        return (foot->getWorldCoordinates(foot->endEffectorPoints[0].pos) + foot->getWorldCoordinates(foot->endEffectorPoints[1].pos)) / 2.0;
    }

    P3D getRightToePosition() {
        pRB foot = getRBByName("right_ankle_roll_link");
        return (foot->getWorldCoordinates(foot->endEffectorPoints[2].pos) + foot->getWorldCoordinates(foot->endEffectorPoints[3].pos)) / 2.0;
    }

    P3D getLeftToePosition() {
        pRB foot = getRBByName("left_ankle_roll_link");
        return (foot->getWorldCoordinates(foot->endEffectorPoints[2].pos) + foot->getWorldCoordinates(foot->endEffectorPoints[3].pos)) / 2.0;
    }

    Quaternion getLeftFootOrientation() {
        pRB foot = getRBByName("left_ankle_roll_link");
        return foot->orientation;
    }

    Quaternion getRightFootOrientation() {
        pRB foot = getRBByName("right_ankle_roll_link");
        return foot->orientation;
    }

    int getLeftAnkleRollIndex() {
        return idx(getJointByName("left_ankle_roll_joint"));
    }

    int getRightAnkleRollIndex() {
        return idx(getJointByName("right_ankle_roll_joint"));
    }

};

typedef std::shared_ptr<H1> pH1;
typedef RobotState H1State;


