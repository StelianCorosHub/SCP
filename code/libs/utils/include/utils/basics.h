#pragma once

// set flag to initialize every eigen matrix with zeros
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#include <float.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#ifndef PI
#define PI 3.1415926535897932
#endif

#ifndef MAX
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#endif

/**
 * This macro checks to see if the value of x is zero within epsilon:
 * -epsilon<x<epsilon
 */
#define IS_ZERO(x) (fabs(x) < 1e-10)

#define IS_EQUAL(x, y) IS_ZERO(((x) - (y)))

#define RAD(x) (((x)*PI) / 180.0)
#define DEG(x) (((x)*180) / PI)

#define SQR(x) ((x) * (x))

#define IS_INFINITE(x) ((x) <= -DBL_MAX || (x) >= DBL_MAX || (x) == INFINITY || (x) == -INFINITY)
#define IS_FINITE(x) ((x) > -DBL_MAX && (x) < DBL_MAX && (x) != INFINITY && (x) != -INFINITY)

#define SGN(x) (((x) < 0) ? (-1) : (1))

// some typedefs
typedef unsigned int uint;
#define Array std::vector

typedef Eigen::Matrix2d Matrix2x2;
typedef Eigen::Matrix3d Matrix3x3;
typedef Eigen::Matrix4d Matrix4x4;
typedef Eigen::Vector2d Vector2d;
typedef Eigen::Vector3d Vector3d;
typedef Eigen::Vector4d Vector4d;

typedef Eigen::VectorXd Vector;
typedef Eigen::MatrixXd Matrix;
typedef Eigen::AngleAxisd AngleAxis;
typedef Eigen::SparseMatrix<double> SparseMatrix;
//sparse matrix triplet: stores value at row i and column j
typedef Eigen::Triplet<double> SMTriplet;

inline void resize(SparseMatrix &sm, int rows, int cols) {
    if (sm.rows() != rows || sm.cols() != cols) sm.resize(rows, cols);
    sm.setZero();
}

inline void resize(Vector &v, int n) {
    if (v.size() != n) v.resize(n);
    v.setZero();
}

inline void resize(Matrix &m, int rows, int cols) {
    if (m.rows() != rows || m.cols() != cols) m.resize(rows, cols);
    m.setZero();
}

inline void clamp(double *v, double min, double max) {
    if (*v < min) *v = min;
    if (*v > max) *v = max;
}

inline void clamp(double &v, double min, double max) {
    if (v < min) v = min;
    if (v > max) v = max;
}

inline void clamp(int &v, int min, int max) {
    if (v < min) v = min;
    if (v > max) v = max;
}

/**
 * if v < min, this method returns 0. If v > max, it returns 1. For everything
 * else it returns an interpolated value;
 */
inline double mapTo01Range(double v, double min, double max) {
    double t = v;
    if (fabs(min - max) < 1e-10) return 1;
    clamp(&t, min, max);
    t = (t - min) / (max - min);
    return t;
}

inline double safeACOS(double val) {
    if (val < -1) return PI;
    if (val > 1) return 0;
    return acos(val);
}

inline double safeASIN(double val) {
    clamp(&val, -1, 1);
    return asin(val);
}

inline double getRandomNumberIn01Range() {
    return static_cast<double>(random()) / RAND_MAX;
}

inline double getRandomNumberInRange(double min, double max) {
    return min + getRandomNumberIn01Range() * (max - min);
}

/**
 * Draw a number from a gaussian distribution. To get a number with a certain
 * mean and variance, take r = mean + sigma*getRandomGaussian()
 */
inline double getRandomNumberFromGaussianDistribution(double m = 0, double s = 1) {
    double x1 = 0, x2 = 0, rquad = 0;
    do {
        x1 = 2.0 * getRandomNumberInRange(0, 1) - 1.0;
        x2 = 2.0 * getRandomNumberInRange(0, 1) - 1.0;
        rquad = x1 * x1 + x2 * x2;
    } while (rquad >= 1 || rquad <= 0);

    double fac = sqrt(-2 * log(rquad) / rquad);

    return m + s * fac * x2;
}

inline bool isNaN(double x) { return (x != x); }

