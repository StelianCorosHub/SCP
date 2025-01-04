
#pragma once

#include <utils/utils.h>
#include <utils/mathUtils.h>

#include <limits>

/**
 * This class is used to represent splines - piecewise polynomial functions
 */
template <class T, class dT> class SplineT {
protected:
    // A caching variable to optimize searching for knots
    mutable int lastIndex = 0;

protected:

//   Splines are defined through control points, which are defined as pairs (t, v), t := knot, or position of the control point, v is its value.
    Array<double> knots;
    Array<T> values;
    Array<dT> tangents;

    /**
     * This method returns the index of the first control point whose knot value is larger
     * than the parameter value t. If no such index exists (i.e. t is larger
     * than any of the values stored), then values.size() is returned.
     */
    inline int getIndexOfFirstControlPointAfter(double t) const {
        int size = (int)knots.size();
        if (size == 0) return 0;

        if (t < knots[(lastIndex + size - 1) % size])
            lastIndex = 0;

        for (int i = 0; i < size; i++) {
            int index = (i + lastIndex) % size;
            if (t < knots[index]) {
                lastIndex = index;
                return index;
            }
        }
        return size;
    }

public:
    SplineT(void) { }

    ~SplineT(void) { }

    /**
     * This method returns a piece-wise constant value of the trajectory at time
     * t (step function)
     */
    inline T eval_piecewiseConstant(double t) const {
        int size = (int)knots.size();
        if (t <= knots[0])
            return values[0];

        if (t >= knots[size - 1])
            return values[size - 1];

        int index = getIndexOfFirstControlPointAfter(t);

        // now figure out where t falls in the correct interval
        t = (t - knots[index - 1]) / (knots[index] - knots[index - 1]);

        if (t < 0.5)
            return values[index - 1];
        else
            return values[index];
    }

    /**
     * This method performs linear interpolation to evaluate the trajectory at
     * the point t
     */
    inline T eval_linear(double t) const {
        int size = (int)knots.size();
        if (size == 0) return zero_init(T());
        if (t <= knots[0])
            return values[0];

        if (t >= knots[size - 1])
            return values[size - 1];

        int index = getIndexOfFirstControlPointAfter(t);

        // now figure out where t falls in the correct interval
        t = (t - knots[index - 1]) / (knots[index] - knots[index - 1]);

        return values[index - 1] * (1-t) + values[index] * t;
    }

    inline T getValRightBefore(double t) const {
        int size = (int)knots.size();
        if (size == 0) return zero_init(T());
        if (t <= knots[0])
            return values[0];
        if (t >= knots[size - 1])
            return values[size - 1];
        int index = getIndexOfFirstControlPointAfter(t);
        return values[index - 1];
    }

    inline T getValRightAfter(double t) const {
        int size = (int)knots.size();
        if (size == 0) return zero_init(T());
        if (t <= knots[0])
            return values[0];
        if (t >= knots[size - 1])
            return values[size - 1];
        int index = getIndexOfFirstControlPointAfter(t);
        return values[index];
    }

    inline double getInterpValue(double t) const {
        int size = (int)knots.size();
        if (size == 0) return 0;
        if (t <= knots[0])
            return 0;

        if (t >= knots[size - 1])
            return 1;

        int index = getIndexOfFirstControlPointAfter(t);

        // now figure out where t falls in the correct interval
        return (t - knots[index - 1]) / (knots[index] - knots[index - 1]);
    }

    /**
     * Evaluate using catmull rom interpolation
     */
    T eval_catmullRom(double t) const {
        int size = (int)knots.size();
        if (knots.size() == 0) return zero_init(T());
        if (t <= knots[0]) return values[0];
        if (t >= knots[size - 1]) return values[size - 1];
        int index = getIndexOfFirstControlPointAfter(t);

        T p1 = values[index - 1];
        T p2 = values[index];

        dT m1 = estimateTangent(index - 1) * (knots[index] - knots[index - 1]);
        dT m2 = estimateTangent(index) * (knots[index] - knots[index - 1]);

        // now that we found the interval, get a value that indicates how far we
        // are along it
        t = (t - knots[index - 1]) / (knots[index] - knots[index - 1]);

        double t2 = t * t, t3 = t * t * t;

        // and now perform the interpolation using the four hermite basis
        // functions from wikipedia
        return p1 * (2 * t3 - 3 * t2 + 1) + m1 * (t3 - 2 * t2 + t) +
               p2 * (-2 * t3 + 3 * t2) + m2 * (t3 - t2);
    }

    /**
 * Evaluate using catmull rom interpolation
 */
    T eval_Hermite(double t) const {
        if (tangents.size() != values.size())
            return eval_catmullRom(t);
        int size = (int)knots.size();
        if (knots.size() == 0) return zero_init(T());
        if (t <= knots[0]) return values[0];
        if (t >= knots[size - 1]) return values[size - 1];
        int index = getIndexOfFirstControlPointAfter(t);

        T p1 = values[index - 1];
        T p2 = values[index];

        dT m1 = tangents[index - 1] * (knots[index] - knots[index - 1]);
        dT m2 = tangents[index] * (knots[index] - knots[index - 1]);

        // now that we found the interval, get a value that indicates how far we
        // are along it
        t = (t - knots[index - 1]) / (knots[index] - knots[index - 1]);

        double t2 = t * t, t3 = t * t * t;

        // and now perform the interpolation using the four hermite basis
        // functions from wikipedia
        return p1 * (2 * t3 - 3 * t2 + 1) + m1 * (t3 - 2 * t2 + t) +
               p2 * (-2 * t3 + 3 * t2) + m2 * (t3 - t2);
    }

    /**
     * Returns the number of knots in this trajectory
     */
    int getNumberOfControlPoints() const { return (int)knots.size(); }

    /**
     * This method is used to insert a new knot in the current trajectory
     */
    void addControlPoint(double t, const T& val) {
        // first we need to know where to insert it, based on the t-values
        int index = getIndexOfFirstControlPointAfter(t);

        knots.insert(knots.begin() + index, t);
        values.insert(values.begin() + index, val);
    }

    /**
     * This method is used to insert a new knot in the current trajectory
     */
    void addControlPoint(double t, const T& val, const dT& tangent) {
        bool useTangents = knots.size() == tangents.size();
        // first we need to know where to insert it, based on the t-values
        int index = getIndexOfFirstControlPointAfter(t);

        knots.insert(knots.begin() + index, t);
        values.insert(values.begin() + index, val);
        if (useTangents)
            tangents.insert(tangents.begin() + index, tangent);
    }

    T getControlPointValue(int idx){
        return values[idx];
    }

    void setControlPointValue(int idx, const T& val) {
        values[idx] = val;
    }

    double getControlPointKnotPosition(int idx) const {
        return knots[idx];
    }

    /**
     * This method is used to remove a knot from the current trajectory.
     * It is assumed that i is within the correct range.
     */
    void removeControlPoint(int i) {
        knots.erase(knots.begin() + i);
        values.erase(values.begin() + i);
    }

    void clear(){
        knots.clear();
        values.clear();
        tangents.clear();
        lastIndex = 0;
    }

    const Array<T>& getValues() const {
        return values;
    }

private:
    dT estimateTangent(int index) const {
        int size = (int)knots.size();
        if (size < 2) return dT();

        if (index == 0)
            return d_(values[0], values[1]) / (knots[1] - knots[0]);

        if (index == size - 1)
            return d_(values[size - 2], values[size - 1]) / (knots[size - 1] - knots[size - 2]);

        dT tBefore = d_(values[index - 1], values[index]) / (knots[index] - knots[index - 1]);
        dT tAfter = d_(values[index], values[index + 1]) / (knots[index + 1] - knots[index]);

        return (tBefore + tAfter) / 2.0;
    }

    double d_(double v1, double v2) const{
        return v2 - v1;
    }

    V3D d_(const P3D& p1, const P3D& p2) const{
        return V3D(p1, p2);
    }

    double zero_init(double val) const{
        return 0;
    }

    P3D zero_init(const P3D& p) const{
        return P3D(0,0,0);
    }

};

typedef SplineT<double, double> Spline;
typedef SplineT<P3D, V3D> Spline3D;



/* old implementation below, in case debugging is needed... */
/*
    template <class T>
    class GenericTrajectory {
    protected:
        // A caching variable to optimize searching for knots
        mutable int lastIndex;

    protected:
        Array<double> tValues;
        Array<T> values;

        inline int getFirstLargerIndex(double t) const {
            int size = (int)tValues.size();
            if (size == 0) return 0;

            if (t < tValues[(lastIndex + size - 1) % size]) lastIndex = 0;
            for (int i = 0; i < size; i++) {
                int index = (i + lastIndex) % size;
                if (t < tValues[index]) {
                    lastIndex = index;
                    return index;
                }
            }
            return size;
        }

    public:
        GenericTrajectory(void) { lastIndex = 0; }

        GenericTrajectory(const GenericTrajectory<T> &other) {
            lastIndex = 0;
            copy(other);
        }

        ~GenericTrajectory(void) { clear(); }


        inline T evaluate_piecewise_constant(double t) const {
            T res;
            evaluate_piecewise_constant(t, res);
            return res;
        }


        inline void evaluate_piecewise_constant(double t, T &res) const {
            int size = (int)tValues.size();
            if (t <= tValues[0]) {
                res = values[0];
                return;
            }
            if (t >= tValues[size - 1]) {
                res = values[size - 1];
                return;
            }
            int index = getFirstLargerIndex(t);

            // now figure out where t falls in the correct interval
            t = (t - tValues[index - 1]) / (tValues[index] - tValues[index - 1]);
            if (t < 0.5)
                res = values[index - 1];
            else
                res = values[index];
        }


        inline T evaluate_linear(double t) const {
            T res;
            evaluate_linear(t, res);
            return res;
        }

        inline void evaluate_linear(double t, T &res) const {
            int size = (int)tValues.size();
            if (t <= tValues[0]) {
                res = values[0];
                return;
            }
            if (t >= tValues[size - 1]) {
                res = values[size - 1];
                return;
            }
            int index = getFirstLargerIndex(t);

            // now linearly interpolate between index-1 and index
            t = (t - tValues[index - 1]) / (tValues[index] - tValues[index - 1]);
            _interp(values[index - 1], values[index], 1 - t, t, res);
        }

        void closest_time_linear(const double &pose, double &t,
                                 std::pair<double, double> timeFrame) const {

            int size = (int)this->values.size();
            boundToRange(timeFrame.first, this->tValues[0],
                         this->tValues[size - 1]);
            boundToRange(timeFrame.second, this->tValues[0],
                         this->tValues[size - 1]);
            int startIndex = this->getFirstLargerIndex(timeFrame.first) - 1;
            int endIndex = this->getFirstLargerIndex(timeFrame.second);
            clamp(startIndex, 0, size - 1);
            clamp(endIndex, 0, size - 1);


            //Get closest index
            int index_min = 0;
            double dist_min = 100000.0;

            for (int i = startIndex; i <= endIndex; i++) {
                double dist_tmp = _distance(this->values[i], pose);
                if (dist_tmp < dist_min) {
                    dist_min = dist_tmp;
                    index_min = i;
                }
            }
            clamp(index_min, 0, size - 1);

            //Get index before or after
            int index_before = index_min - 1;
            clamp(index_before, 0, size - 1);
            int index_after = index_min + 1;
            clamp(index_after, 0, size - 1);
            double dist_before = _distance(this->values[index_before], pose);
            double dist_after = _distance(this->values[index_after], pose);

            int index_other = index_min;
            double dist_other = dist_min;
            if (dist_before < dist_after) {
                index_other = index_before;
                dist_other = dist_before;
            } else {
                index_other = index_after;
                dist_other = dist_after;
            }

            //Linear interpolation of time
            double normalize = dist_min + dist_other;
            if (index_min == index_other || normalize < 1e-6) {
                t = this->tValues[index_min];
                return;
            }
            t = (1.0 - (dist_other / normalize)) * this->tValues[index_other] +
                (1.0 - (dist_min / normalize)) * this->tValues[index_min];
        }

        T getSlopeEstimateAtKnot(int index,
                                 bool equalEndpointSlopes = false) const {
            if (getKnotCount() < 2) return T();

            if (index == 0 || index == getKnotCount() - 1) {
                T startSlope = static_cast<T>((values[1] - values[0]) /
                                              (tValues[1] - tValues[0]));
                T endSlope = static_cast<T>(
                        (values[getKnotCount() - 1] - values[getKnotCount() - 2]) /
                        (tValues[getKnotCount() - 1] - tValues[getKnotCount() - 2]));

                if (equalEndpointSlopes) return (startSlope + endSlope) / 2.0;

                if (index == 0) return startSlope;
                return endSlope;
            }

            T slopeBefore = static_cast<T>((values[index] - values[index - 1]) /
                                           (tValues[index] - tValues[index - 1]));
            T slopeAfter = static_cast<T>((values[index + 1] - values[index]) /
                                          (tValues[index + 1] - tValues[index]));

            return (slopeBefore + slopeAfter) / 2.0;
        }

        double length() const {
            int size = tValues.size();
            double t1 = tValues[size - 1];
            double t0 = tValues[0];

            double numSteps = 10;
            double dt = (t1 - t0) / (numSteps - 1.0);
            double t = t0;
            double l = 0;
            for (int i = 0; i < numSteps; ++i) {
                l += evaluate_catmull_rom(t);
                t += dt;
            }
            return l;
        }

        T evaluate_catmull_rom(double t, bool equalEndpointSlopes = false) const {
            int size = (int)tValues.size();
            if (tValues.size() == 0) return T(0);
            if (t <= tValues[0]) return values[0];
            if (t >= tValues[size - 1]) return values[size - 1];
            int index = getFirstLargerIndex(t);

            // now that we found the interval, get a value that indicates how far we
            // are along it
            t = (t - tValues[index - 1]) / (tValues[index] - tValues[index - 1]);

            // approximate the derivatives at the two ends

            T p1 = values[index - 1];
            T p2 = values[index];

            T m1 = getSlopeEstimateAtKnot(index - 1, equalEndpointSlopes) *
                   (tValues[index] - tValues[index - 1]);
            T m2 = getSlopeEstimateAtKnot(index, equalEndpointSlopes) *
                   (tValues[index] - tValues[index - 1]);

            double t2, t3;
            t2 = t * t;
            t3 = t2 * t;

            // and now perform the interpolation using the four hermite basis
            // functions from wikipedia
            return p1 * (2 * t3 - 3 * t2 + 1) + m1 * (t3 - 2 * t2 + t) +
                   p2 * (-2 * t3 + 3 * t2) + m2 * (t3 - t2);
        }

        T getKnotValue(int i) const { return values[i]; }

        double getKnotPosition(int i) const { return tValues[i]; }

        void setKnotValue(int i, const T &val) { values[i] = val; }

        void setKnotPosition(int i, double pos) {
            if (i - 1 >= 0 && tValues[i - 1] >= pos) return;
            if ((uint)(i + 1) < tValues.size() - 1 && tValues[i + 1] <= pos) return;
            tValues[i] = pos;
        }

        double getMinPosition() const {
            if (tValues.empty()) return std::numeric_limits<double>::infinity();
            return tValues.front();
        }

        double getMaxPosition() const {
            if (tValues.empty()) return -std::numeric_limits<double>::infinity();
            return tValues.back();
        }


        int getKnotCount() const { return (int)tValues.size(); }


        void addKnot(double t, T val) {
            // first we need to know where to insert it, based on the t-values
            int index = getFirstLargerIndex(t);

            tValues.insert(tValues.begin() + index, t);
            values.insert(values.begin() + index, val);
        }


        void removeKnot(int i) {
            tValues.erase(tValues.begin() + i);
            values.erase(values.begin() + i);
        }

        void clear() {
            tValues.clear();
            values.clear();
        }

        void copy(const GenericTrajectory<T> &other) {
            tValues.clear();
            values.clear();
            int size = other.getKnotCount();

            tValues.reserve(size);
            values.reserve(size);
            for (int i = 0; i < size; ++i) {
                tValues.push_back(other.tValues[i]);
                values.push_back(other.values[i]);
            }
        }

    private:
        inline static void _interp(double p1, double p2, double t1, double t2,
                                   double &res) {
            res = p1 * t1 + p2 * t2;
        }

        inline static void _interp(const V3D &v1, const V3D &v2, double t1,
                                   double t2, V3D &res) {
            res[0] = v1[0] * t1 + v2[0] * t2;
            res[1] = v1[1] * t1 + v2[1] * t2;
            res[2] = v1[2] * t1 + v2[2] * t2;
        }

        inline static void _interp(const P3D &p1, const P3D &p2, double t1,
                                   double t2, P3D &res) {
            res[0] = p1[0] * t1 + p2[0] * t2;
            res[1] = p1[1] * t1 + p2[1] * t2;
            res[2] = p1[2] * t1 + p2[2] * t2;
        }

        inline static void _interp(const Eigen::VectorXd &p1,
                                   const Eigen::VectorXd &p2, double t1, double t2,
                                   Eigen::VectorXd &res) {
            res = p1 * t1 + p2 * t2;
        }

        inline static double _distance(const Eigen::VectorXd &p1,
                                       const Eigen::VectorXd &p2) {
            return (p1 - p2).norm();
        }

        inline static double _distance(const double &p1, const double &p2) {
            return fabs(p1 - p2);
        }
    };

    typedef GenericTrajectory<double> Trajectory1D;
    typedef GenericTrajectory<V3D> Trajectory3D;
*/

