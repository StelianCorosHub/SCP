#pragma once

#include <utils/mathUtils.h>
#include <utils/utils.h>
#include <utils/physicsUtils.h>
#include <RBSim/RBUtils.h>

class ArticulatedFigure;

/**
 * This class acts as a container for the state of an articulated figure
 **/
class AFPose {
public:
    Quaternion rootQ;
    P3D rootPos;

    //joint angles
    Array<double> q;

public:
    AFPose() {}

    AFPose(int nJoints){
        q.resize(nJoints, 0);
    }

    AFPose(const ArticulatedFigure& af);
    AFPose(ArticulatedFigure* af);
    AFPose(std::shared_ptr<ArticulatedFigure> pAF);

    ~AFPose() {}

    inline AFPose(const AFPose &other) {
        rootQ = other.rootQ;
        rootPos = other.rootPos;

        q = other.q;
    }

    inline double getHeadingAngle() const {
        return rootQ.getHeadingAngle(PHYS::up);
    }

    //setting the heading angle will rotate the entire state, including linear and
    //angular velocities, around the world up axis. The root position will remain unaffected.
    void setHeadingAngle(double headingAngle) {
        // get the current heading
        Quaternion oldHeading = rootQ.getHeading(PHYS::up);
        // this is now the heading that we want
        Quaternion newHeading = Quaternion(headingAngle, PHYS::up) * oldHeading.inverse();
        // add this component to the root.
        rootQ = newHeading * rootQ;
    }

    void writeToFile(const char *fName) {
        if (fName == nullptr)
            throwError("cannot write to a file whose name is nullptr!");

        FILE *fp = fopen(fName, "w");

        if (fp == nullptr)
            throwError("cannot open file \'%s\' for writing...", fName);

        fprintf(fp, "%lf %lf %lf\n", rootPos.x, rootPos.y, rootPos.z);
        fprintf(fp, "%lf %lf %lf %lf\n", rootQ.w(), rootQ.x(), rootQ.y(), rootQ.z());

        fprintf(fp, "%d\n", (int)q.size());

        for (int i = 0; i < q.size(); i++)
            fprintf(fp, "%lf\n", q[i]);

        fclose(fp);
    }

    void readFromFile(const char *fName) {
        if (fName == nullptr)
            throwError("cannot read a file whose name is nullptr!");

        FILE *fp = fopen(fName, "r");

        if (fp == nullptr)
            throwError("cannot open file \'%s\' for reading...", fName);

        fscanf(fp, "%lf %lf %lf", &rootPos.x, &rootPos.y, &rootPos.z);
        double t1 = 1, t2 = 0, t3 = 0, t4 = 0;
        fscanf(fp, "%lf %lf %lf %lf", &t1, &t2, &t3, &t4);
        rootQ.w() = t1; rootQ.x() = t2; rootQ.y() = t3; rootQ.z() = t4;

        int n = 0;
        fscanf(fp, "%d", &n);
        q.resize(n);

//        for (int i = 0; i < q.size(); i++) {
//            fscanf(fp, "%lf %lf %lf %lf", &t1, &t2, &t3, &t4);
//            Quaternion rotq(t1, t2, t3, t4);
//            q[i] = rotq.getRotationAngle();
//        }

        for (int i = 0; i < q.size(); i++) {
            fscanf(fp, "%lf", &t1);
            q[i] = t1;
        }

        fclose(fp);
    }

    bool operator == (const AFPose& other) {
        if (q.size() != other.q.size())
            return false;

        double errMax = 1e-7;

        if (V3D(rootPos, other.rootPos).norm() > errMax)
            return false;
        if ((rootQ.inverse() * other.rootQ).vec().norm() > errMax)
            return false;

        for (uint i = 0; i<q.size(); i++) {
            if (fabs(q[i] - other.q[i]) > errMax)
                return false;
        }
        return true;
    }
};

/**
 * This class acts as a container for the state of an articulated figure
 **/
class AFState {
public:
    Quaternion rootQ;
    P3D rootPos;
    V3D rootVel;
    V3D rootAngVel;

    //joint angles
    Array<double> q;
    //rate of change of joint angles
    Array<double> qDot;

public:
    AFState() {}

    AFState(int nJoints) {
        q.resize(nJoints, 0);
        qDot.resize(nJoints, 0);
    }

    AFState(const AFPose &afp) {
        rootQ = afp.rootQ;
        rootPos = afp.rootPos;
        rootVel = V3D();
        rootAngVel = V3D();

        q = afp.q;
        qDot = q;
        for (int i = 0; i < q.size(); i++) qDot[i] = 0;
    }

    AFState(const AFPose &afp, const AFPose &afp2_dt, double dt) {
        rootQ = afp.rootQ;
        rootPos = afp.rootPos;
        rootVel = V3D(rootPos, afp2_dt.rootPos) / dt;
        rootAngVel = estimateAngularVelocity(rootQ, afp2_dt.rootQ, dt);

        q = afp.q; qDot = q;
        for (int i = 0; i < qDot.size(); i++)
            qDot[i] = (afp2_dt.q[i] - q[i]) / dt;
    }

    AFState(const ArticulatedFigure& af);
    AFState(ArticulatedFigure* af);
    AFState(std::shared_ptr<ArticulatedFigure> pAF);

    ~AFState() {}

    inline AFState(const AFState &other) {
        rootQ = other.rootQ;
        rootPos = other.rootPos;
        rootVel = other.rootVel;
        rootAngVel = other.rootAngVel;

        q = other.q;
        qDot = other.qDot;
    }

    void estimateQDotBasedOnDifferenceTo(const AFState& other, double dt);

    inline double getHeadingAngle() const {
        return rootQ.getHeadingAngle(PHYS::up);
    }

    //TODO: test this - is it the case that if you integrate the angular velocity forward in time, recompute heading,
    // do finite differences, then we get the same result?
    inline double getRateOfChangeOfHeading() const {
        return rootAngVel.dot(PHYS::up);
    }

    //setting the heading angle will rotate the entire state, including linear and
    //angular velocities, around the world up axis. The root position will remain unaffected.
    void setHeadingAngle(double headingAngle) {
        // get the current heading
        Quaternion oldHeading = rootQ.getHeading(PHYS::up);
        // this is now the heading that we want
        Quaternion newHeading = Quaternion(headingAngle, PHYS::up) * oldHeading.inverse();
        // add this component to the root.
        rootQ = newHeading * rootQ;
        rootVel = newHeading * rootVel;
        rootAngVel = newHeading * rootAngVel;
    }

    void writeToFile(const char *fName) {
        if (fName == nullptr)
            throwError("cannot write to a file whose name is nullptr!");

        FILE *fp = fopen(fName, "w");

        if (fp == nullptr)
            throwError("cannot open file \'%s\' for writing...", fName);

        fprintf(fp, "%lf %lf %lf\n", rootPos.x, rootPos.y, rootPos.z);
        fprintf(fp, "%lf %lf %lf %lf\n", rootQ.w(), rootQ.x(), rootQ.y(), rootQ.z());
        fprintf(fp, "%lf %lf %lf\n", rootVel[0], rootVel[1], rootVel[2]);
        fprintf(fp, "%lf %lf %lf\n\n", rootAngVel[0], rootAngVel[1], rootAngVel[2]);

        fprintf(fp, "%d\n", (int)q.size());

        for (int i = 0; i < q.size(); i++)
            fprintf(fp, "%lf %lf\n", q[i], qDot[i]);

        fclose(fp);
    }

    void readFromFile(const char *fName){
        if (fName == nullptr)
            throwError("cannot read a file whose name is nullptr!");

        FILE *fp = fopen(fName, "r");

        if (fp == nullptr)
            throwError("cannot open file \'%s\' for reading...", fName);

        fscanf(fp, "%lf %lf %lf", &rootPos.x, &rootPos.y, &rootPos.z);
        double t1 = 1, t2 = 0, t3 = 0, t4 = 0;
        fscanf(fp, "%lf %lf %lf %lf", &t1, &t2, &t3, &t4);
        rootQ.w() = t1; rootQ.x() = t2; rootQ.y() = t3; rootQ.z() = t4;
        fscanf(fp, "%lf %lf %lf", &rootVel[0], &rootVel[1], &rootVel[2]);
        fscanf(fp, "%lf %lf %lf", &rootAngVel[0], &rootAngVel[1], &rootAngVel[2]);

        int n = 0;
        fscanf(fp, "%d", &n);
        q.resize(n);
        qDot.resize(n);

//        for (int i = 0; i < jq.size(); i++) {
//            fscanf(fp, "%lf %lf %lf %lf", &t1, &t2, &t3, &t4);
//            Quaternion qrot(t1, t2, t3, t4);
//            jq[i] = q.getRotationAngle(); jqDot[i] = 0;
//            fscanf(fp, "%lf %lf %lf", &t1, &t2, &t3);
//            ... set qdot here...
//        }

        for (int i = 0; i < q.size(); i++){
            fscanf(fp, "%lf %lf", &t1, &t2);
            q[i] = t1; qDot[i] = t2;
        }

        fclose(fp);
    }

    bool operator == (const AFState& other) {
        if (q.size() != other.q.size())
            return false;

        double errMax = 1e-7;

        if (V3D(rootPos, other.rootPos).norm() > errMax)
            return false;
        if ((rootVel - other.rootVel).norm() > errMax)
            return false;
        if ((rootAngVel - other.rootAngVel).norm() > errMax)
            return false;

        if ((rootQ.inverse() * other.rootQ).vec().norm() > errMax)
            return false;
        for (uint i = 0; i<q.size(); i++) {
            if (fabs(q[i] - other.q[i]) > errMax)
                return false;
            if (fabs(qDot[i] - other.qDot[i]) > errMax)
                return false;
        }
        return true;
    }
};

