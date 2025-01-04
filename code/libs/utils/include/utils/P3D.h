#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

class P3D {
public:
    double x = 0;
    double y = 0;
    double z = 0;

    explicit P3D(double x, double y, double z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    //equivalent to *this = P3D() + v;
    explicit P3D(const Eigen::Vector3d& v){
        this->x = v.x();
        this->y = v.y();
        this->z = v.z();
    }

    P3D() { this->x = this->y = this->z = 0; }

    explicit P3D(const double *data) {
        this->x = data[0];
        this->y = data[1];
        this->z = data[2];
    }

    P3D operator+(const P3D &p) const { return P3D(x + p.x, y + p.y, z + p.z); }

    P3D operator-(const P3D &p) const { return P3D(x - p.x, y - p.y, z - p.z); }

    P3D &operator += (const P3D &p) {
        x += p.x;
        y += p.y;
        z += p.z;
        return *this;
    }

    P3D &operator -= (const P3D &p) {
        x -= p.x;
        y -= p.y;
        z -= p.z;
        return *this;
    }

    P3D operator * (double v) const { return P3D(x * v, y * v, z * v); }

    P3D &operator *= (double v) {
        x *= v;
        y *= v;
        z *= v;
        return *this;
    }

    P3D operator / (double v) { return P3D(x / v, y / v, z / v); }

    P3D &operator /= (double v) {
        x /= v;
        y /= v;
        z /= v;
        return *this;
    }

    double &operator [] (int idx) {
        if (idx == 0) return x;
        if (idx == 1) return y;
        return z;
    }

    const double &operator [] (int idx) const {
        if (idx == 0) return x;
        if (idx == 1) return y;
        return z;
    }

    bool isAtInfinity(){
        return IS_INFINITE(x) || IS_INFINITE(y) || IS_INFINITE(z);
    }
};

