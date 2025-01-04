#pragma once

#include <Eigen/Dense>

#include <utils/P3D.h>

class V3D : public Eigen::Vector3d {
public:
    explicit V3D(double x, double y, double z) {
        (*this)[0] = x;
        (*this)[1] = y;
        (*this)[2] = z;
    }

    explicit V3D() {
        (*this)[0] = 0;
        (*this)[1] = 0;
        (*this)[2] = 0;
    }

    /**
     * this constructor will get called by default every time Vector3d
     * operations are involved
     */
    V3D(const Eigen::Vector3d &v) { Eigen::Vector3d::operator=(v); }

    /**
     * vector that points from p1 to p2
     */
    V3D(const P3D &p1, const P3D &p2) {
        (*this)[0] = p2.x - p1.x;
        (*this)[1] = p2.y - p1.y;
        (*this)[2] = p2.z - p1.z;
    }

    /**
     * vector that points from origin to p
     */
    explicit V3D(const P3D &p) {
        (*this)[0] = p.x;
        (*this)[1] = p.y;
        (*this)[2] = p.z;
    }

    V3D &operator = (const Eigen::Vector3d &v) {
        Eigen::Vector3d::operator=(v);
        return *this;
    }

    V3D operator + (const V3D &v) const { return V3D(Eigen::Vector3d::operator+(v)); }

    V3D operator + (const Eigen::Vector3d &v) const {
        return V3D(Eigen::Vector3d::operator+(v));
    }

    V3D operator - () const { return (V3D)Eigen::Vector3d::operator-(); }

    V3D operator - (const V3D &v) const { return (V3D)Eigen::Vector3d::operator-(v); }

    V3D operator - (const Eigen::Vector3d &v) const {
        return (V3D)Eigen::Vector3d::operator-(v);
    }

    V3D operator * (double val) const { return V3D(Eigen::Vector3d::operator*(val)); }

    V3D operator / (double val) const { return V3D(Eigen::Vector3d::operator/(val)); }

    V3D cross(const V3D &other) const { return V3D(Eigen::Vector3d::cross(other)); }

    double getComponentAlong(const V3D &other) { return Eigen::Vector3d::dot(other); };

    V3D unit() const {
        if (this->norm() < 10e-20) return V3D(1, 0, 0);
        return *this / this->norm();
    }

    //returns a unit vector (one of the inifinite possibilities) that is orthogonal to *this
    V3D getOrthonormalVector() const {
        if (this->norm() < 10e-20) return V3D(1, 0, 0);
        V3D guess = V3D(1,0,0);
        //check if *this and our guess are almost aligned. If so, guess another one that is
        // orthogonal to our first guess. They can't both be aligned with *this.
        if (this->dot(guess) < this->norm() * 0.01)
            guess = V3D(0,0,1);
        return this->cross(guess).unit();
    }

    double getComponentAlong(const V3D &other) const {
        return Eigen::Vector3d::dot(other);
    };

    void setComponentAlong(const V3D &other, double val) {
        double oldVal = getComponentAlong(other);
        *this += other * (val - oldVal);
    };

    double angleWith(const Eigen::Vector3d& other) const {
        // U.V = |U|*|V|*cos(angle)
        // therefore angle = inverse cos (U.V/(|U|*|V|))
        if (other.norm() < 1e-10 || this->norm() < 1e-10)
            return 0;
        return safeACOS(this->dot(other) / (other.norm() * this->norm()));
    }

    double signedAngleWith(const Eigen::Vector3d& other, const Eigen::Vector3d& positiveAxis) const {
        double angle = this->angleWith(other);
        if (this->cross(other).dot(positiveAxis) < 0)
            angle *= -1;
        return angle;
    }

    bool isZero(){
        return norm() <= 1e-10;
    }
};

// return val * v, a V3D
inline V3D operator * (double val, const V3D &v) { return v * val; }

// return p + v, a P3D
inline P3D operator + (const P3D &p, const Eigen::Vector3d &v) {
    return P3D(p.x + v[0], p.y + v[1], p.z + v[2]);
}

// return p + v, a P3D
inline P3D &operator += (P3D &p, const Eigen::Vector3d &v) {
    p.x += v[0];
    p.y += v[1];
    p.z += v[2];
    return p;
}

//return p - v, a P3D
inline P3D operator - (const P3D &p, const Eigen::Vector3d &v) {
    return P3D(p.x - v[0], p.y - v[1], p.z - v[2]);
}
