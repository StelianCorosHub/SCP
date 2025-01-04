#include "RBSim/RBSRenderer.h"

#include <gui/renderer.h>

#define NORMAL_COLOR V3D(0.8, 0.8, 0.8)
#define HIGHLIGHTED_COLOR V3D(1.0, 0.5, 0.5)
#define COM_COLOR V3D(0.5, 1.0, 0.5)

#define SKELETON_RADIUS 0.01
#define COLLISION_PRIMITIVE_COLOR V3D(0.75, 0.0, 0.0)

#define EE_DRAW_COLOR V3D(0.0, 0.0, 1.0)
#define CONTACT_EE_DRAW_COLOR V3D(0.0, 1.0, 1.0)

#define JOINT_AXIS_COLOR V3D(1.0, 0.5, 0.5)
#define JOINT_COLOR V3D(0.5, 0.5, 0.5)

void RBSRenderer::drawSkeleton(const pRB& rb, const pRBJ pJoint, const Array<pRBJ> cJoints, double alpha) {
    // draw capsules that define the "skeleton" of this body: parent joints TO
    // features TO child joints and end effectors

    V3D drawColor = rb->selected ? HIGHLIGHTED_COLOR : NORMAL_COLOR;

    drawSphere(rb->position, SKELETON_RADIUS * 1.5, COM_COLOR, alpha);

    P3D centerP = P3D(0,0,0);
    int n = cJoints.size();
    if (pJoint.get()) n += 1;
    if (n >= 2) {
        if (pJoint.get())
            centerP = pJoint->cJPos;
        for (uint i = 0; i < cJoints.size(); i++)
            centerP += cJoints[i]->pJPos;
        centerP /= n;
    }

    drawCapsule(rb->position, rb->getWorldCoordinates(centerP), SKELETON_RADIUS * 0.2,
                drawColor, alpha);

    if (pJoint != NULL) {
        P3D startPos = rb->getWorldCoordinates(pJoint->cJPos);
        P3D endPos = rb->getWorldCoordinates(centerP);
        drawCapsule(startPos, endPos, SKELETON_RADIUS, drawColor, alpha);
        drawSphere(startPos, SKELETON_RADIUS * 1.3, drawColor, alpha);
    }

    for (uint i = 0; i < cJoints.size(); i++) {
        P3D startPos = rb->getWorldCoordinates(centerP);
        P3D endPos = rb->getWorldCoordinates(cJoints[i]->pJPos);
        drawCapsule(startPos, endPos, SKELETON_RADIUS, drawColor, alpha);
        drawSphere(endPos, SKELETON_RADIUS * 1.3, drawColor, alpha);
    }

    if (cJoints.size() == 0) {
        for (uint i = 0; i < rb->endEffectorPoints.size(); i++) {
            P3D startPos = rb->getWorldCoordinates(centerP);
            P3D endPos = rb->getWorldCoordinates(rb->endEffectorPoints[i].pos);
            drawCapsule(startPos, endPos, SKELETON_RADIUS, drawColor, alpha);
        }
    }
}

void RBSRenderer::draw3DModels(const pRB& rb, double alpha) {
    for (uint i = 0; i < rb->RBMeshes.size(); i++) {
        RigidTransform meshTransform(rb->orientation, rb->position);
        meshTransform *= rb->RBMeshes[i]->transform;

        rb->RBMeshes[i]->model->position = meshTransform.T;
        rb->RBMeshes[i]->model->orientation = meshTransform.R;

        if (rb->selected)
            rb->RBMeshes[i]->model->draw(HIGHLIGHTED_COLOR, alpha);
        else
            rb->RBMeshes[i]->model->draw(rb->RBMeshes[i]->color, alpha);
    }
}

void RBSRenderer::drawCollisionPrimitives(const pRB& rb) {
    for (uint i = 0; i < rb->collisionShapes.size(); i++) {
        rb->collisionShapes[i]->draw(COLLISION_PRIMITIVE_COLOR, 0.5);
    }
}

void RBSRenderer::drawMOI(const Matrix3x3& moiLocal, double mass, const P3D& pos, const Quaternion& rot, double alpha, bool wireFrame){
    //first of all, draw the COM as a green sphere
//    drawSphere(pos, 0.02, V3D(0,1,0));

    //triple check that the moi is a symmetric matrix
    assert(IS_ZERO(moiLocal(1,0)-moiLocal(0,1)) &&
                   IS_ZERO(moiLocal(2,0)-moiLocal(0,2))  &&
                   IS_ZERO(moiLocal(1,2)-moiLocal(2,1)) );

    Eigen::EigenSolver<Matrix3x3> eigenvalueSolver(moiLocal, true);

    Eigen::Vector3cd principleMomentsOfInertia = eigenvalueSolver.eigenvalues();

    assert(IS_ZERO(principleMomentsOfInertia[0].imag()) &&
           IS_ZERO(principleMomentsOfInertia[1].imag()) &&
           IS_ZERO(principleMomentsOfInertia[1].imag()));

    Matrix3x3 V = eigenvalueSolver.eigenvectors().real();

    double Ixx = principleMomentsOfInertia[0].real();  // = m(y2 + z2)/12
    double Iyy = principleMomentsOfInertia[1].real();  // = m(z2 + x2)/12
    double Izz = principleMomentsOfInertia[2].real();  // = m(y2 + x2)/12

    double x = sqrt((Iyy + Izz - Ixx) * 6 / mass);
    double y = sqrt((Izz + Ixx - Iyy) * 6 / mass);
    double z = sqrt((Ixx + Iyy - Izz) * 6 / mass);

    P3D pmin(-x / 2, -y / 2, -z / 2), pmax(x / 2, y / 2, z / 2);

    if (V.determinant() < 0.0) {
        V(0, 2) *= -1;
        V(1, 2) *= -1;
        V(2, 2) *= -1;
    }

    //eigenvectors of eigenvalues with multiplicity > 1 need not be orthogonal... so do a bit of a projection here...
    V.col(1) = V.col(1) - V.col(0).dot(V.col(1)) * V.col(0); V.col(1).normalize();

    V.col(2) = V.col(2) - V.col(0).dot(V.col(2)) * V.col(0); V.col(2).normalize();
    V.col(2) = V.col(2) - V.col(1).dot(V.col(2)) * V.col(1); V.col(2).normalize();

    assert(IS_ZERO(abs(V.determinant() - 1.0)) &&
           "Rotation matrices have a determinant which is equal to 1.0!");

    Quaternion q(V.real());

    if (wireFrame == false)
        drawCuboid(pos, rot * q, V3D(x, y, z),
                   V3D(0.7, 0.7, 0.7), alpha);
    else
        drawWireFrameCuboid(pos, rot * q, V3D(x, y, z),
                            V3D(1.0, 0.7, 0.7));
}

void RBSRenderer::drawMOI(const pRB& rb, double alpha) {
    drawMOI(rb->MOILocal(), rb->Mass(), rb->position, rb->orientation, alpha, false);
}

void RBSRenderer::drawEndEffectors(const pRB& rb) {
    for (uint i = 0; i < rb->endEffectorPoints.size(); i++) {
        auto &ee = rb->endEffectorPoints[i];
        P3D pos = rb->getWorldCoordinates(ee.pos);
        V3D color = ee.inContact ? CONTACT_EE_DRAW_COLOR
                                 : EE_DRAW_COLOR;
        drawSphere(pos, ee.radius, color);
    }
}

void RBSRenderer::drawJoint(const pRBJ& j) {
    drawArrow3d(j->getWorldPosition(),V3D(j->parent->getWorldCoordinates(j->rotAxis) * 0.1), 0.01,
                JOINT_AXIS_COLOR);

    drawArrow3d(j->getWorldPosition(),V3D(j->child->getWorldCoordinates(j->rotAxis) * 0.1), 0.01,
                JOINT_AXIS_COLOR);

    drawSphere(j->parent->getWorldCoordinates(j->pJPos), SKELETON_RADIUS * 1.3, JOINT_COLOR);
    drawSphere(j->child->getWorldCoordinates(j->cJPos), SKELETON_RADIUS * 1.3, JOINT_COLOR);
}


