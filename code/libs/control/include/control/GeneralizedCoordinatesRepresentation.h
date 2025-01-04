#pragma once

#include <unordered_map>
#include <RBSim/RB.h>
#include <RBSim/ArticulatedFigure.h>
#include <utils/mathUtils.h>

/**
 * This class implements a reduced coordinate representation for articulated rigid body systems
 * (e.g. robots). Their configuration will be represented via the generalized coordinate vectors q and qDot
 */
class GeneralizedCoordinatesRepresentation {
private:
    //--- Reduced-coordinates state representation
    Vector q, qDot;

    //For each DOF, store its rotation/translation axis, expressed in local coordinates
    Array<V3D> qAxisLocal;
    //index of the parent DOF for each q
    Array<int> qParentIdx;
    Array<V3D> qOffsetToParent;
    Array<V3D> qOffsetToChildLinkCOM;

public:
    /** the constructor */
    GeneralizedCoordinatesRepresentation(ArticulatedFigure *af);

    GeneralizedCoordinatesRepresentation() { }

    /** the destructor */
    virtual ~GeneralizedCoordinatesRepresentation(void) {}

    /**
     * updates q and qDot given the state of the articulated figure
     */
    void setGeneralizedCoordinateValuesFromAFState(const AFState& afState);

    // updates q and qDot given the pose of the articulated figure
    void setGeneralizedCoordinateValuesFromAFPose(const AFPose& afPose);

    /**
     * sets the state of the articulated figure based on current q / qDot values
     */
    AFState getAFState();

    /**
     * sets the pose of the articulated figure based on current q values
     */
    AFPose getAFPose();

    int getQSize() { return q.size(); }

    Vector getQ() { return q; }
    Vector getQDot() { return qDot; }

    void setQ(const Vector& q) { assert(q.size() == this->q.size()); this->q = q;}
    void setQDot(const Vector& qDot) { assert(qDot.size() == this->qDot.size()); this->qDot = qDot;}

    // returns the world coordinates for point p, which is specified in the local
    // coordinates of the child link (relative to its COM) of qIdx: p(q)
    P3D p_of_q(const P3D &pLocal, int qIdx);

    // returns the world-coordinates velocity for point pLocal, which is specified
    // in the local coordinates of the child link of qIdx (relative to its COM)
    V3D pDot(const P3D &pLocal, int qIdx);

    // returns the world coordinates for vector v, which is specified in the local
    // coordinates of the child link of qIdx: v(q)
    V3D v_of_q(const V3D &vLocal, int qIdx);

    // returns the global orientation for child link of qIdx: R(q)
    Quaternion R_of_q(int qIdx);

    // returns the world-coordinates angular velocity of the child link of qIdx
    V3D RDot(int qIdx);



    /**
     * returns the axis corresponding to the indexed generalized coordinate,
     * expressed in local coordinates.
     */
    V3D getQAxis(int qIndex) const;

    /**
     * returns the translation or rotation axis for a specific dof q...
     */
    V3D getQAxis_world(int qIndex);

    /**
    * In the tree-like hierarchy of joints/dofs, this method returns the parent index of the
    * dof corresponding to qIdx
    */
    int getParentQIdxOf(int qIdx);

    void projectVectorOnGeneralizedCoordsAxes(const V3D &vector, const V3D &a,
                                              const V3D &b, const V3D &c,
                                              double &aVal, double &bVal,
                                              double &cVal);


    /**
     * computes the jacobian dp/dq that tells you how the world coordinates of point p
     * change with q. p is expressed in the local coordinates of the child link of qIdx
     */
    void compute_dpdq(const P3D &pLocal, int qIdx, Matrix &dpdq);

    /**
     * computes the jacobian dp/dq that tells you how the world coordinates of vector v
     * change with q. v is expressed in the local coordinates of the child link of qIdx
     */
     void compute_dvdq(const V3D &vLocal, int qIdx, Matrix &dvdq);

    // computes the jacobian dR/dq which relates changes in the
    // orientation (represented in axis-angle form) of the child link of qIdx to changes in q
    void compute_dRdq(int qIdx, Matrix &dRdq);


    /**
     * estimates the jacobian dp/dq using finite differences
     */
    void estimate_linear_jacobian(const P3D &pLocal, int qIdx, Matrix &dpdq);

    /**
     * estimates the jacobian dv/dq using finite differences
     */
    void estimate_linear_jacobian(const V3D &vLocal, int qIdx, Matrix &dvdq);

    bool test_linear_jacobian(const P3D &pLocal, int qIdx);
    bool test_linear_jacobian(const V3D &vLocal, int qIdx);


    void estimate_angular_jacobian(int qIdx, Matrix &dRdq);
    bool test_angular_jacobian(int qIdx);

    /**
     * computes the matrix that tells you how the jacobian dp/dq changes with
     * respect to q_i. Returns true if it contains non-zero elements, false
     * otherwise
     */
    bool compute_ddpdq_dqi(const P3D &pLocal, int qIdx, Matrix &ddpdq_dqi, int q_i);

    /**
     * computes the matrix that tells you how the jacobian dv/dq changes with
     * respect to q_i. Returns true if it contains non-zero elements, false
     * otherwise
     */
    bool compute_ddvdq_dqi(const V3D &vLocal, int qIdx, Matrix &ddvdq_dqi, int q_i);

    /**
     * estimates the change of dp/dq with respect to q_i
     */
    void estimate_ddpdq_dqi(const P3D &pLocal, int qIdx, Matrix &ddpdq_dqi, int q_i);

    /**
     * estimates the change of dv/dq with respect to q_i
     */
    void estimate_ddvdq_dqi(const V3D &vLocal, int qIdx, Matrix &ddvdq_dqi, int q_i);

    bool test_linear_jacobian_derivatives(const V3D &vLocal, int qIdx);
    bool test_linear_jacobian_derivatives(const P3D &vLocal, int qIdx);

    // computes d(Jw)/dqi
    bool compute_ddRdq_dqi(int qIdx, Matrix &ddRdqdqi, int q_i);

    // estimates the change of angular jacobian with respect to q_i using finite
    // differences
    void estimate_ddRdq_dqi(int qIdx, Matrix &ddRdq_dqi, int q_i);

    bool test_angular_jacobian_derivatives(int qIdx);

};

void testGeneralizedCoordinateRepresentation(Robot *robot);

