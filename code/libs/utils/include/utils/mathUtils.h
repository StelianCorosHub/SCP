#pragma once

#include <utils/basics.h>
#include <utils/P3D.h>
#include <utils/V3D.h>

#define UP V3D(0,1,0)

#include <utils/Quaternion.h>
#include <utils/RigidTransform.h>

/**
 * this returns the angle between u and v, result is in the range 0..PI
 */
inline double angleBetween(const V3D &u, const V3D &v) {
    return u.angleWith(v);
}

/**
 * tells us the rotation angle between u and v, given a rotation axis n
 */
inline double angleBetween(const V3D &u, const V3D &v, const V3D &n) {
    double a = angleBetween(u, v);
    if (u.cross(v).dot(n) < 0) a = 2 * PI - a;
    return a;
}

/**
 * returns skew symmetric matrix corresponding to v x
 */
inline Matrix3x3 getSkewSymmetricMatrix(const V3D &v) {
    Matrix3x3 result;
    result << 0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0;
    return result;
}

inline Matrix3x3 getCrossProductMatrix(const V3D &v) {
    return getSkewSymmetricMatrix(v);
}

inline void getOrthonormalVectors(const V3D& n, V3D& p, V3D& q) {

    if (fabs(n[2]) > M_SQRT1_2) {
        // choose p in y-z plane
        double a = n[1]*n[1] + n[2]*n[2];
        double k = 1 / sqrt(a);
        p[0] = 0;
        p[1] = -n[2]*k;
        p[2] = n[1]*k;
        // set q = n x p
        q[0] = a*k;
        q[1] = -n[0]*p[2];
        q[2] = n[0]*p[1];
    }
    else {
        // choose p in x-y plane
        double a = n[0]*n[0] + n[1]*n[1];
        double k = 1 / sqrt(a);
        p[0] = -n[1]*k;
        p[1] = n[0]*k;
        p[2] = 0;
        // set q = n x p
        q[0] = -n[2]*p[1];
        q[1] = n[2]*p[0];
        q[2] = a*k;
    }
}