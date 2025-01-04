#pragma once

#include <utils/mathUtils.h>

class RigidTransform {
public:
    Quaternion R;
    P3D T = P3D(0, 0, 0);

public:
    RigidTransform(const Quaternion &_R = Quaternion(),
                        const P3D &_T = P3D())
            : R(_R), T(_T) {}

    ~RigidTransform() {}

    P3D transform(const P3D &p) { return T + R * V3D(p); }

    V3D transform(const V3D &v) { return R * v; }

    RigidTransform inverse() {
        RigidTransform trans;
        trans.R = R.inverse();
        trans.T = P3D() - (trans.R * V3D(P3D(), T));
        return trans;
    }

    RigidTransform operator*(const RigidTransform &other) {
        RigidTransform trans;

        trans.R = R * other.R;
        trans.T = T + R * V3D(other.T);
        return trans;
    }

    RigidTransform &operator*=(const RigidTransform &other) {
        RigidTransform trans;

        trans.R = R * other.R;
        trans.T = T + R * V3D(other.T);
        *this = trans;
        return *this;
    }
};

