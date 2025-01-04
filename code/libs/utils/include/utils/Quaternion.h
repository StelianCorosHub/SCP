#pragma once

#include <Eigen/Dense>
#include <utils/basics.h>

class Quaternion : public Eigen::Quaternion<double> {
public:

    //default constructor which creates an identity quaternion
    Quaternion() {
        (*this).x() = 0;
        (*this).y() = 0;
        (*this).z() = 0;
        (*this).w() = 1;
    }

    //creates a quaternion. Angle is in radians, axis is a unit
    Quaternion(double angle, const Eigen::Vector3d& axis){
        *this = Eigen::Quaternion<double>(Eigen::AngleAxisd(angle, axis));
    }

    //creates a quaternion from an angle axis rotation representation
    explicit Quaternion(const Eigen::AngleAxisd& aa){
        *this = Eigen::Quaternion<double>(aa);
    }

    //creates a quaternion from a rotation matrix
    explicit Quaternion(const Eigen::Matrix3d& R){
        *this = Eigen::Quaternion<double>(R);
    }

    //creates the smallest quaternion that aligns vectors a and b (i.e. a = *this * b, b = this->inverse() * a)
    explicit Quaternion(const V3D& a, const V3D& b) {
        V3D rotAxis = a.cross(b);
        if (rotAxis.norm() > 1e-12)
            *this = Eigen::Quaternion<double>(Eigen::AngleAxisd(-safeACOS(a.dot(b) / (a.norm() * b.norm())), rotAxis.unit()));
        else {
            //a and b are either aligned (an identity quaternion will do), or they are
            // rotated relative to each other by 180 about any of the many axes that are
            // orthogonal to both a and b
            if (a.dot(b) > 0)
                *this = Quaternion();
            else
                *this = Eigen::Quaternion<double>(Eigen::AngleAxisd(PI, a.getOrthonormalVector()));
        }
    }

    Quaternion(double w, double x, double y, double z) {
        (*this).x() = x;
        (*this).y() = y;
        (*this).z() = z;
        (*this).w() = w;
    }

    Quaternion &operator = (const Eigen::Quaternion<double> &q) {
        Eigen::Quaternion<double>::operator=(q);
        return *this;
    }

    /**
     * this constructor will get called by default every time Quaternion
     * operations are involved
     */
    Quaternion(const Eigen::Quaternion<double> &q) { Eigen::Quaternion<double>::operator=(q); }

    double getRotationAngle(const Eigen::Vector3d& v) const {
        double result = 2 * safeACOS(w());
        if (vec().dot(v) < 0) result = -result;
        if (result > PI) result -= 2 * PI;
        if (result < -PI) result += 2 * PI;
        return result;
    }

    double getRotationAngle() const {
        return getRotationAngle(vec().normalized());
    }

/**
 * decompose the quaternion q as: q = R(c, gamma) * R(b, beta) * R(a, alpha).
 * Unknowns are: alpha, beta, gamma.
 * Assumptions: a, b and c are unit vectors; b is orthogonal to a and c; a is
 * either orthogonal to c, or is aligned with it.
 */
    void computeEulerAngles(const V3D &a, const V3D &b, const V3D &c,
                            double &alpha, double &beta, double &gamma) const {
        // the idea here is that the a axis only gets rotated about b and c, which
        // are assumed to be orthogonal to each other. Based on this info, we can
        // first compute the angles beta and gamma
        assert(IS_ZERO(a.dot(b)) && IS_ZERO(b.dot(c)));
        assert(IS_ZERO(a.norm() - 1) && IS_ZERO(b.norm() - 1) &&
               IS_ZERO(c.norm() - 1));

        V3D aRot = *this * a;

        if (IS_ZERO(a.dot(c))) {
            bool circular = a.cross(b).dot(c) > 0;
            // the three axes form an orthonormal basis (i.e. Tait-Bryan)...
            // singularity around beta = -PI/2 or PI/2
            if (circular) {
                beta = -safeASIN(aRot.dot(c));
                gamma = atan2(aRot.dot(b), aRot.dot(a));
            } else {
                beta = safeASIN(aRot.dot(c));
                gamma = atan2(-aRot.dot(b), aRot.dot(a));
            }
        } else if (IS_ZERO(a.dot(c) - 1)) {
            // these are "proper" euler axes, where the first and the last one are
            // the same... singularity around beta = 0 or PI
            V3D lastAxis = a.cross(b);
            beta = safeACOS(aRot.dot(a));
            gamma = atan2(aRot.dot(b), -aRot.dot(lastAxis));
        } else {
            // dunno what this is.... freak out...
            alpha = beta = gamma = 0;
            assert(false);
        }

        Quaternion qLeft = Quaternion(-beta, b) * Quaternion(-gamma, c) * *this;
        alpha = qLeft.getRotationAngle(a);

        double tmpVal = (Quaternion(gamma, c) * Quaternion(beta, b) * Quaternion(alpha, a) * this->inverse()).vec().norm();
        if (tmpVal >= 1e-3){
            double tmpVal2 = this->norm();
            assert(tmpVal < 1e-3);
        }

//        assert((Quaternion(gamma, c) * Quaternion(beta, b) * Quaternion(alpha, a) * this->inverse()).vec().norm() < 1e-3);
    }

    /**
     * This method decomposes the current rotation quaternion (pQc) into a
     * heading / yaw / rotation about the UP axis component and everything
     * else (e.g. mix of roll and pitch). It returns the first
     * component.
     *
     * Think of the current quaternion as representing the relative orientation
     * between two coordinate frames p and c (i.e. q rotates vectors from a
     * local / child frame c into a global / parent frame p: pQc). With u = UP
     * specified in the global coordinate frame p, this method decomposes
     * pQc such that:
     *
     *      pQc = pQt * tQc, where pQt is purely a rotation about axis u.
     *
     * Note that in the intermediate coordinate frame t, u has the same
     * coordinates as in the global frame of reference p, and tQc is the
     * rotation that aligns the coordinates of the vector u from c to those
     * from coordinate frames t and p.
     */

    Quaternion getHeading(const V3D& u = UP) const {
        // we need to compute u in c's coordinates
        V3D u_c = this->inverse() * u;

        // now we compute the smallest rotation that aligns u in the parent
        // and child coordinate frames
        Quaternion tQc(u, u_c);

        Quaternion pQt = *this * tQc.inverse();

        return pQt;
    }

    // returns the angle of the heading part of the rotation, expressed in radians
    double getHeadingAngle(const V3D& u = UP) const {
        return getHeading(u).getRotationAngle(u);
    }
};
