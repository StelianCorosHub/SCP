#pragma once

#include <utils/logger.h>

#include <utils/BVHReader.h>

#include "guiMath.h"
#include "model.h"

inline void drawSphere(const P3D &p, double r,
                       const V3D &color = V3D(1, 0, 0), double alpha = 1.0) {
    static Model sphere = Model(SCP_DATA_FOLDER "/meshes/sphere.obj");

    sphere.position = p;
    sphere.scale = V3D(2 * r, 2 * r, 2 * r);
    sphere.draw(color, alpha);
}

/**
 * p is center of ellipsode
 * orientation is orientation of ellipsoid
 * dims is radius along each axis (in ellipsoid's frame)
 */
inline void drawEllipsoid(const P3D &p, const Quaternion &orientation,
                          const V3D &dims,
                          const V3D &color = V3D(1, 0, 0), double alpha = 1.0) {
    static Model sphere = Model(SCP_DATA_FOLDER "/meshes/sphere.obj");

    sphere.position = p;
    sphere.orientation = orientation;
    sphere.scale = V3D(2 * dims.x(), 2 * dims.y(), 2 * dims.z());
    sphere.draw(color, alpha);
}

inline void drawCuboid(const P3D &p, const Quaternion &orientation,
                       const V3D &dims,
                       const V3D &color = V3D(1, 0, 0), double alpha = 1.0) {
    static Model cube = Model(SCP_DATA_FOLDER "/meshes/cube.obj");

    cube.position = p;
    cube.orientation = orientation;
    cube.scale = dims;
    cube.draw(color, alpha);
}

inline void drawWireFrameCuboid(const P3D &p, const Quaternion &orientation,
                                const V3D &dims,
                                const V3D &color = V3D(1, 0, 0),
                                double alpha = 1.0) {
    static Model cube = Model(SCP_DATA_FOLDER "/meshes/cube_frame.obj");

    cube.position = p;
    cube.orientation = orientation;
    cube.scale = dims;
    cube.draw(color, alpha);
}

inline void drawCylinder(const P3D &startPosition, const P3D &endPosition,
                         const double &radius,
                         const V3D &color = V3D(1, 0, 0), double alpha = 1.0) {
    static Model cylinder = Model(SCP_DATA_FOLDER "/meshes/cylinder.obj");

    V3D dir = V3D(startPosition, endPosition);
    double s = dir.norm();
    if (s < 10e-10) return;
    V3D a = dir.normalized();
    V3D b = V3D(0, 0, 1);
    V3D v = (b.cross(a)).normalized();
    // can happen that dir and the z-vector are directly aligned, in which case
    // the angle will be 0 or PI, but we need a valid axis...
    if (v.norm() < 0.5) v = V3D(1, 0, 0);
    float angle = acos(b.dot(a) / (b.norm() * a.norm()));

    cylinder.position = startPosition;
    cylinder.scale = V3D(radius, radius, s);
    cylinder.orientation = Quaternion(AngleAxis(angle, v));
    cylinder.draw(color, alpha);
}

inline void drawCylinder(const P3D &startPosition, const V3D &direction,
                         const double &radius,
                         const V3D &color = V3D(1, 0, 0), double alpha = 1.0) {
    drawCylinder(startPosition, startPosition + direction, radius,
                 color, alpha);
}

inline void drawCone(const P3D &origin, const V3D &direction,
                     const double &radius,
                     const V3D &color = V3D(1, 0, 0), double alpha = 1.0) {
    static Model cone = Model(SCP_DATA_FOLDER "/meshes/cone.obj");

    // Transformation & Scale
    double s = direction.norm();
    if (s < 10e-10) return;
    V3D a = direction.normalized();
    V3D b = V3D(0, 1, 0);
    V3D v = (b.cross(a)).normalized();
    // can happen that dir and the y-vector are directly aligned, in which case
    // the angle will be 0 or PI, but we need a valid axis...
    if (v.norm() < 0.5) v = V3D(1, 0, 0);
    float angle = acos(b.dot(a) / (b.norm() * a.norm()));

    cone.position = origin;
    cone.scale = V3D(1e-3 * radius, 1e-3 * s, 1e-3 * radius);
    cone.orientation = Quaternion(AngleAxis(angle, v));
    cone.draw(color, alpha);
}

inline void drawArrow3d(const P3D &origin, const V3D &direction,
                        const double &radius,
                        const V3D &color = V3D(1, 0, 0), double alpha = 1.0) {
    V3D dir = direction;
    for (uint i = 0; i < 3; i++)
        if (fabs(dir[i]) < 1e-6) dir[i] = 1e-6;

    double coneRadius = 1.5 * radius;
    V3D coneDir = dir / dir.norm() * coneRadius * 1.5;
    P3D cyl_endPos = origin + dir - coneDir;

    drawCylinder(origin, cyl_endPos, radius, color, alpha);
    drawCone(cyl_endPos, coneDir, coneRadius, color, alpha);
}

inline void drawCapsule(const P3D &sP, const P3D &eP, const double &radius,
                        const V3D &color = V3D(1, 0, 0),
                        double alpha = 1.0) {
    drawCylinder(sP, eP, 2*radius, color, alpha);
    drawSphere(sP, radius, color, alpha);
    drawSphere(eP, radius, color, alpha);
}

inline void drawCapsule(const P3D &startPosition, const V3D &direction,
                        const double &radius,
                        const V3D &color = V3D(1, 0, 0), double alpha = 1.0) {
    drawCapsule(startPosition, P3D(startPosition + direction), radius,
                color, alpha);
}

inline void drawBVHPose(const MocapSkeletonPose& sp, double radius = 0.1, double alpha = 1.0, int selectedBIdx = -1) {

    for (int i = 0; i < sp.getBoneCount(); i++) {
        V3D col = V3D(0.6, 0.6, 0.6);
        if (i == selectedBIdx)
            col = V3D(1.0, 0.6, 0.6);
        if (V3D(sp.getBoneStartPos(i), sp.getBoneEndPos(i)).norm() > 2 * radius)
            drawCapsule(sp.getBoneStartPos(i), sp.getBoneEndPos(i), radius * 0.8, col, alpha);
    }

    for (int i = 0; i < sp.getSFCount(); i++) {
        if (sp.getSFType(i) == SF_JOINT)
            drawSphere(sp.getSFPos(i), radius * 0.99, V3D(1,0,0), alpha);
        if (sp.getSFType(i) == SF_END_EFFECTOR)
            drawSphere(sp.getSFPos(i), radius, V3D(0,0,1), alpha);
        if (sp.getSFType(i) == SF_ROOT){
//            drawArrow3d(sp.getSFPos(i), sp.getRootOrientation() * (V3D(1,0,0) * 5 * radius), radius * 0.1, V3D(1,0,0));
//            drawArrow3d(sp.getSFPos(i), sp.getRootOrientation() * (V3D(0,1,0) * 5 * radius), radius * 0.1, V3D(0,1,0));
//            drawArrow3d(sp.getSFPos(i), sp.getRootOrientation() * (V3D(0,0,1) * 5 * radius), radius * 0.1, V3D(0,0,1));
            drawCuboid(sp.getSFPos(i), sp.getRootOrientation(), V3D(3 * radius, 3 * radius, 3 * radius), V3D(0,0,0), 0.25 * alpha);
        }
    }

}

inline void drawCircle(const P3D &center, const V3D &normal,
                       const double &radius,
                       const V3D &color = V3D(1, 0, 0), double alpha = 1.0) {
    // Compute start and end position
    double height = 1e-3;
    drawCylinder(center, center + (normal.normalized() * height), radius,
                 color, alpha);
}

inline void drawRectangle(const P3D &p, const V3D &normal, const double &angle,
                          const Vector2d &dims,
                          const V3D &color = V3D(1, 0, 0), double alpha = 1.0) {
    V3D dim;
    dim[0] = dims[0];
    dim[1] = 1e-3;
    dim[2] = dims[1];
    V3D a = normal.normalized();
    V3D b = V3D(0, 1, 0);
    V3D v = (b.cross(a)).normalized();
    float theta = acos(b.dot(a) / (b.norm() * a.norm()));

    drawCuboid(p,
               Quaternion(AngleAxis(theta, v)) *
                   Quaternion(AngleAxis(angle, V3D(0, 1, 0))),
               dim, color, alpha);
}

// we could go from the vector "from" the long way, or the short way. The Vector
// up will tell us which one is meant
inline void drawSector(const P3D &p, const V3D &from, const V3D &to,
                       const V3D &up,
                       const V3D &color = V3D(1, 0, 0), double alpha = 1.0) {
    static Model sector(SCP_DATA_FOLDER "/meshes/sector.obj");

    double radius = (from.norm() + to.norm()) / 2.0;
    double angleMax = angleBetween(from, to, up);

    drawArrow3d(p, from, radius / 40.0, color, alpha = 1.0);
    drawArrow3d(p, to, radius / 40.0, color, alpha = 1.0);

    // Transformation & Scale
    V3D a = from.cross(to).normalized();
    if (a.dot(up) < 0) a *= -1;

    sector.position = p;
    sector.scale = V3D(1e-3 * radius, 1e-3 * radius, 1e-3 * radius);
    Matrix3x3 rot;

    // the x axis will map to the vector "from" -- in the mesh coordinate frame,
    // x corresponds to the axis the arc is on
    rot.col(0) = from.normalized();
    // the y axis will map to the normal...
    rot.col(1) = a;
    // and the z axis will be the cross product of the two...
    rot.col(2) = from.cross(a).normalized();

    Quaternion qDefault(rot);

    sector.orientation = qDefault;

    // Draw
    for (double ang = 0; ang < angleMax; ang += PI / 180.0) {
        sector.orientation =
            qDefault * Quaternion(AngleAxis(ang + PI / 2.0, V3D(0, 1, 0)));
        sector.draw(color, alpha);
    }
}


