#pragma once

#include <gui/model.h>
#include <gui/shader.h>

#include "RBSim/RBJoint.h"

class RBSRenderer {
public:

    static void drawSkeleton(const pRB& rb, const pRBJ pJoint, const Array<pRBJ> cJoints, double alpha = 1.0);

    static void draw3DModels(const pRB& rb, double alpha = 1.0);

    static void drawCollisionPrimitives(const pRB& rb);

    static void drawMOI(const Matrix3x3& moi, double mass, const P3D& pos, const Quaternion& q, double alpha, bool wireFrame);

    static void drawMOI(const pRB& rb, double alpha = 0.5);

    static void drawEndEffectors(const pRB& rb);

    static void drawJoint(const pRBJ& j);
};
