#include <control/GeneralizedCoordinatesRepresentation.h>
#include <RBSim/RBJoint.h>
#include <utils/utils.h>

// TODO: for this implementation, we have chosen not to precompute any
// quantities (e.g. world coords of joints or joint axes, position/orientations
// of RBs, etc), trading speed in favor of a conceptually clean implementation.
// If ever the need will arise, things like RB orientations, positions and
// rotation axes of joints, etc can be precomputed as soon as the q's get set,
// which should make things faster, but at the cost of additional management
// of data.


//TODO: do rethink the q for rotation of the root - get best parameterization, not x y z...

GeneralizedCoordinatesRepresentation::GeneralizedCoordinatesRepresentation(ArticulatedFigure *af) {
    resize(q, 6 + (int)af->joints.size());
    resize(qDot, 6 + (int)af->joints.size());

    //first 6 dofs correspond to the 3 translation + 3 rotation DOFs of the root
    qAxisLocal.push_back(V3D(1,0,0)); qParentIdx.push_back(-1); qOffsetToParent.push_back(V3D(0,0,0)); qOffsetToChildLinkCOM.push_back(V3D(0,0,0));
    qAxisLocal.push_back(V3D(0,1,0)); qParentIdx.push_back(0); qOffsetToParent.push_back(V3D(0,0,0)); qOffsetToChildLinkCOM.push_back(V3D(0,0,0));
    qAxisLocal.push_back(V3D(0,0,1)); qParentIdx.push_back(1); qOffsetToParent.push_back(V3D(0,0,0)); qOffsetToChildLinkCOM.push_back(V3D(0,0,0));

    //use a roll pitch yaw representation here, though this could change easily when needed...
    //y - yaw
    qAxisLocal.push_back(V3D(0, 1, 0)); qParentIdx.push_back(2); qOffsetToParent.push_back(V3D(0,0,0)); qOffsetToChildLinkCOM.push_back(V3D(0,0,0));
    // x - pitch
    qAxisLocal.push_back(V3D(1, 0, 0)); qParentIdx.push_back(3); qOffsetToParent.push_back(V3D(0,0,0)); qOffsetToChildLinkCOM.push_back(V3D(0,0,0));
    // z - roll
    qAxisLocal.push_back(V3D(0, 0, 1)); qParentIdx.push_back(4); qOffsetToParent.push_back(V3D(0,0,0)); qOffsetToChildLinkCOM.push_back(V3D(0,0,0));

    //after the first 6dofs, we will use the exact same indexing for the joints as in the representation of the articulated figure...
    int jIdx = 0;
    for (auto j : af->joints) {
        qAxisLocal.push_back(j->rotAxis);

        if (af->getParentOf(j->parent) == nullptr)
            qOffsetToParent.push_back(V3D(P3D(), j->pJPos));
        else
            qOffsetToParent.push_back(V3D(af->getParentOf(j->parent)->cJPos, P3D()) + V3D(P3D(), j->pJPos));

        qParentIdx.push_back(6 + af->idx(af->getParentOf(j->parent)));
        qOffsetToChildLinkCOM.push_back(V3D(j->cJPos, P3D()));
    }

    setGeneralizedCoordinateValuesFromAFState(AFState(af));
}

// updates q and qDot given the state of the robot
void GeneralizedCoordinatesRepresentation::setGeneralizedCoordinateValuesFromAFState(const AFState& afState) {
    q[0] = afState.rootPos.x;
    q[1] = afState.rootPos.y;
    q[2] = afState.rootPos.z;

    afState.rootQ.computeEulerAngles(getQAxis(5), getQAxis(4),getQAxis(3), q[5], q[4], q[3]);

    projectVectorOnGeneralizedCoordsAxes(afState.rootVel, getQAxis(0), getQAxis(1), getQAxis(2), qDot[0], qDot[1], qDot[2]);

    projectVectorOnGeneralizedCoordsAxes(afState.rootAngVel, getQAxis_world(3), getQAxis_world(4), getQAxis_world(5), qDot[3], qDot[4], qDot[5]);

    // Now go through each joint, and decompose it as appropriate...
    for (uint i = 0; i < afState.q.size(); i++) {
        q[6 + i] = afState.q[i];
        qDot[6 + i] = afState.qDot[i];
    }
}

// updates q and qDot given the pose of the robot
void GeneralizedCoordinatesRepresentation::setGeneralizedCoordinateValuesFromAFPose(const AFPose& afPose) {
    q[0] = afPose.rootPos.x;
    q[1] = afPose.rootPos.y;
    q[2] = afPose.rootPos.z;

    afPose.rootQ.computeEulerAngles(getQAxis(5), getQAxis(4),getQAxis(3), q[5], q[4], q[3]);

    qDot[0] = qDot[1] = qDot[2] = qDot[3] = qDot[4] = qDot[5] = 0;

    // Now go through each joint, and decompose it as appropriate...
    for (uint i = 0; i < afPose.q.size(); i++) {
        q[6 + i] = afPose.q[i];
        qDot[6 + i] = 0;
    }
}

/**
 * sets the state of the articulated figure based on current q / qDot values
 */
AFState GeneralizedCoordinatesRepresentation::getAFState() {
    AFState state(q.size() - 6);

    // set the position, velocity, rotation and angular velocity for the root
    state.rootPos = P3D() + getQAxis(0) * q[0] + getQAxis(1) * q[1] + getQAxis(2) * q[2];
    state.rootVel = getQAxis(0) * qDot[0] + getQAxis(1) * qDot[1] + getQAxis(2) * qDot[2];
    state.rootAngVel =  getQAxis_world(3) * qDot[3] +
                        getQAxis_world(4) * qDot[4] +
                        getQAxis_world(5) * qDot[5];
    state.rootQ = R_of_q(5);

    for (uint i = 0; i < state.q.size(); i++) {
        state.q[i] = q[6 + i];
        state.qDot[i] = qDot[6 + i];
    }

    // and done...
    return state;
}

/**
 * sets the pose of the articulated figure based on current q values
 */
AFPose GeneralizedCoordinatesRepresentation::getAFPose() {
    AFPose pose(q.size() - 6);

    // set the position, velocity, rotation and angular velocity for the root
    pose.rootPos = P3D() + getQAxis(0) * q[0] + getQAxis(1) * q[1] + getQAxis(2) * q[2];
    pose.rootQ = R_of_q(5);

    for (uint i = 0; i < pose.q.size(); i++)
        pose.q[i] = q[6 + i];

    // and done...
    return pose;
}

/**
* In the tree-like hierarchy of joints/dofs, this method returns the parent index of the
* dof corresponding to qIdx
*/
int GeneralizedCoordinatesRepresentation::getParentQIdxOf(int qIdx) {
    return qParentIdx[qIdx];
}

// returns the axis corresponding to the indexed generalized coordinate,
// expressed in local coordinates
V3D GeneralizedCoordinatesRepresentation::getQAxis(int qIdx) const {
    return qAxisLocal[qIdx];
}

// returns the world coordinates for point p, which is specified in the local
// coordinates of the child link (relative to its COM) of qIdx: p(q)
P3D GeneralizedCoordinatesRepresentation::p_of_q(const P3D &pLocal, int qIdx) {
    //this is the offset vector from the location of qIdx to point p, in the local coordinate frame of the child link of qIdx
    V3D vLocal = qOffsetToChildLinkCOM[qIdx] + V3D(P3D(), pLocal);

    // 2 here is the index of the first translational DOF of the root
    while (qIdx > 2) {
        vLocal = qOffsetToParent[qIdx] + Quaternion(q[qIdx], getQAxis(qIdx)) * vLocal;

        //and repeat... compute the local coordinates of the point in the frame of the parent after joint rotation is applied
        qIdx = getParentQIdxOf(qIdx);
    }

    return P3D() + V3D(getQAxis(0) * q[0] + getQAxis(1) * q[1] + getQAxis(2) * q[2]) + vLocal;
}

// computes the jacobian dp/dq that tells you how the world coordinates of p
// change with q. p is expressed in the local coordinates of the child link of qIdx.
void GeneralizedCoordinatesRepresentation::compute_dpdq(const P3D &pLocal, int qIdx, Matrix &dpdq) {
    resize(dpdq, 3, q.size());
    dpdq.fill(0);

    // todo: this bit of code (and other jacobians/higher order derivatives) can
    // probably be made faster if ever they become the bottleneck... seems like
    // it would be easy enough to precompute orientations/joint
    // positions/offsets in a first pass from the bottom up, instead of always
    // recomputing them for each j... after all, we know that dpdq_i is the
    // cross product of the vector from pos of q_i to p (all in world coords)
    // with the rotation axis of q_i (in world coords)...

    int loopIndex = qIdx;
    // 2 here is the index of the first translational DOF of the root
    while (loopIndex > 2) {
        //this is the point in the local coordinate frame of rb
        V3D vLocal = qOffsetToChildLinkCOM[qIdx] + V3D(P3D(), pLocal);

        int qIndex = qIdx;

        // compute the offset from the location of joint qIndex to point p - if
        // these quantitites were known in world coordinates, then it would be
        // clear we wouldn't need the loop
        while (qIndex > loopIndex) {
            //we need to compute the coordinates of pLocal in the frame of the parent DOF
            vLocal = qOffsetToParent[qIndex] + Quaternion(q[qIndex], getQAxis(qIndex)) * vLocal;

            //and repeat... compute the local coordinates of the point in the frame of the parent after joint rotation is applied
            qIndex = getParentQIdxOf(qIndex);
        }

        // this is how the vector computed above changes with the rotation angle
        V3D offset = getQAxis(qIndex).cross(vLocal);

        // and now bring this vector to world coordinates, since the step above
        // was computed in the coordinate frame of joint q. again, if RB
        // orientations were precomputed, we could just read it off now...
        while (qIndex > 2) {
            offset = Quaternion(q[qIndex], getQAxis(qIndex)) * offset;
            qIndex = getParentQIdxOf(qIndex);
        }

        for (int i = 0; i < 3; i++) dpdq(i, loopIndex) = offset[i];

        loopIndex = getParentQIdxOf(loopIndex);
    }

    dpdq(0, 0) = 1;
    dpdq(1, 1) = 1;
    dpdq(2, 2) = 1;
}

// returns the world coordinates for vector v, which is specified in the local
// coordinates of the child link of qIdx: v(q)
V3D GeneralizedCoordinatesRepresentation::v_of_q(const V3D &vLocal, int qIdx) {
    V3D v = vLocal;
    // 2 here is the index of the first translational DOF of the root
    while (qIdx > 2) {
        v = V3D(Quaternion(q[qIdx], getQAxis(qIdx)) * v);
        qIdx = getParentQIdxOf(qIdx);
    }

    return v;
}

// computes the jacobian dv/dq that tells us how the world coordinates of v
// change with q. vLocal is expressed in the local coordinates of the child link of qIdx
void GeneralizedCoordinatesRepresentation::compute_dvdq(const V3D &vLocal, int qIdx, Matrix &dvdq) {
    resize(dvdq, 3, q.size());
    dvdq.fill(0);

    int startIndex = qIdx;

    int loopIndex = startIndex;
    // 2 here is the index of the first translational DOF of the root
    while (loopIndex > 2) {
        V3D theVector(vLocal);
        int qIndex = startIndex;

        while (qIndex > loopIndex) {
            theVector = Quaternion(q[qIndex], getQAxis(qIndex)) * theVector;
            qIndex = getParentQIdxOf(qIndex);
        }

        theVector = getQAxis(qIndex).cross(theVector);

        while (qIndex > 2) {
            theVector = Quaternion(q[qIndex], getQAxis(qIndex)) * theVector;
            qIndex = getParentQIdxOf(qIndex);
        }

        for (int i = 0; i < 3; i++) dvdq(i, loopIndex) = theVector[i];

        loopIndex = getParentQIdxOf(loopIndex);
    }
}

// returns the global orientation for child link of qIdx: R(q)
Quaternion GeneralizedCoordinatesRepresentation::R_of_q(int qIdx) {
    Quaternion qRes = Quaternion::Identity();
    // 2 here is the index of the first translational DOF of the root -- these
    // dofs do not contribute to the orientation of the rigid bodies...
    while (qIdx > 2) {
        qRes = Quaternion(q[qIdx], getQAxis(qIdx)) * qRes;
        qIdx = getParentQIdxOf(qIdx);
    }
    return qRes;
}

// computes the jacobian dR/dq which relates changes in the
// orientation (represented in axis-angle form) of the child link of qIdx to changes in q
void GeneralizedCoordinatesRepresentation::compute_dRdq(int qIdx, Matrix &dRdq) {
    // the orientation of an RB is obtained by rotating by all the rotation axes
    // up the hierarchy... so the jacobian consists of the world-coords axes
    resize(dRdq, 3, (int)q.size());

    std::vector<int> qIndices;

    while (qIdx > 2) {
        qIndices.push_back(qIdx);
        qIdx = getParentQIdxOf(qIdx);
    }

    Quaternion Q = Quaternion::Identity();
    // q rotates the rigid body (grand...child) about its world coordinate
    // axis... these will be the entries of dR/dq...
    while (qIndices.size() > 0) {
        int index = qIndices.back();
        if (qIndices.back() > 2)
            Q = Q * Quaternion(q[qIndices.back()], getQAxis(qIndices.back()));
        // this is just the rotation axis in world coordinates, as this is what
        // we are rotating about to get the orientation of the rigid body if
        // everything else is fixed...
        V3D dRdq_i = Q * getQAxis(index);
        dRdq(0, index) = dRdq_i[0];
        dRdq(1, index) = dRdq_i[1];
        dRdq(2, index) = dRdq_i[2];

        qIndices.pop_back();
    }
}

void GeneralizedCoordinatesRepresentation::estimate_angular_jacobian(int qIdx, Matrix &dRdq) {
    resize(dRdq, 3, (int)q.size());

    for (int i = 0; i < q.size(); i++) {
        double val = q[i];
        double h = 0.001;

        q[i] = val + h;
        Quaternion R_p = R_of_q(qIdx);

        q[i] = val - h;
        Quaternion R_m = R_of_q(qIdx);

        q[i] = val;
        Quaternion rotate = R_p * R_m.inverse();
        Eigen::AngleAxisd aa(rotate);
        V3D axis = aa.axis().normalized();
        double angle = aa.angle();
        axis *= angle;

        V3D dRdq_i = axis / (2 * h);
        dRdq(0, i) = dRdq_i[0];
        dRdq(1, i) = dRdq_i[1];
        dRdq(2, i) = dRdq_i[2];
    }
}

bool GeneralizedCoordinatesRepresentation::test_angular_jacobian(int qIdx) {
    Matrix dRdq_analytic, dRdq_estimated;
    compute_dRdq(qIdx, dRdq_analytic);
    estimate_angular_jacobian(qIdx, dRdq_estimated);

    bool error = false;

    for (int i = 0; i < dRdq_analytic.rows(); i++)
        for (int j = 0; j < dRdq_analytic.cols(); j++) {
            double err = dRdq_analytic(i, j) - dRdq_estimated(i, j);
            if (fabs(err) > 0.0001) {
                std::cout << "Error at: " << i << " " << j
                          << ": analytic: " << dRdq_analytic(i, j)
                          << " estimated " << dRdq_estimated(i, j)
                          << " error: " << err << std::endl;
                error = true;
            }
        }

    return !error;
}

void GeneralizedCoordinatesRepresentation::
    projectVectorOnGeneralizedCoordsAxes(const V3D &vec, const V3D &a,
                                         const V3D &b, const V3D &c,
                                         double &aVal, double &bVal,
                                         double &cVal) {
    // we cannot assume the vectors form an orthogonal basis, so we have to
    // solve a system to compute the unknowns aVal, bVal, cVal. We assume that
    // vec, a, b and c are all specified in the same coordinate frame
    Matrix3x3 m, mInv;
    m(0, 0) = a[0];
    m(0, 1) = b[0];
    m(0, 2) = c[0];
    m(1, 0) = a[1];
    m(1, 1) = b[1];
    m(1, 2) = c[1];
    m(2, 0) = a[2];
    m(2, 1) = b[2];
    m(2, 2) = c[2];

    mInv = m.inverse();

    V3D sol = V3D(mInv * vec);

    aVal = sol[0];
    bVal = sol[1];
    cVal = sol[2];

    // check the solution...
    V3D res = V3D(m * sol) - vec;
    if (res.norm() > 0.0001)
        std::cout << "Failed to properly compute generalized coordinates for "
                     "the input :("
                  << std::endl;
}

// this is a somewhat slow function to use if we must iterate through multiple
// rigid bodies...
V3D GeneralizedCoordinatesRepresentation::getQAxis_world(int qIndex) {
    if (qIndex < 3) return getQAxis(qIndex);
    return R_of_q(qIndex) * getQAxis(qIndex);
}

// estimates the linear jacobian dp/dq using finite differences
void GeneralizedCoordinatesRepresentation::estimate_linear_jacobian(const P3D &pLocal, int qIdx, Matrix &dpdq) {
    resize(dpdq, 3, (int)q.size());

    for (int i = 0; i < q.size(); i++) {
        double val = q[i];
        double h = 0.0001;

        q[i] = val + h;
        P3D p_p = p_of_q(pLocal, qIdx);

        q[i] = val - h;
        P3D p_m = p_of_q(pLocal, qIdx);

        q[i] = val;
        V3D dpdq_i = V3D(p_m, p_p) / (2 * h);
        dpdq(0, i) = dpdq_i[0];
        dpdq(1, i) = dpdq_i[1];
        dpdq(2, i) = dpdq_i[2];
    }
}

// estimates the jacobian dv/dq using finite differences
void GeneralizedCoordinatesRepresentation::estimate_linear_jacobian(const V3D &vLocal, int qIdx, Matrix &dvdq) {
    resize(dvdq, 3, (int)q.size());

    for (int i = 0; i < q.size(); i++) {
        double val = q[i];
        double h = 0.0001;

        q[i] = val + h;
        V3D p_p = v_of_q(vLocal, qIdx);

        q[i] = val - h;
        V3D p_m = v_of_q(vLocal, qIdx);

        q[i] = val;
        V3D dvdq_i = (p_p - p_m) / (2 * h);
        dvdq(0, i) = dvdq_i[0];
        dvdq(1, i) = dvdq_i[1];
        dvdq(2, i) = dvdq_i[2];
    }
}

bool GeneralizedCoordinatesRepresentation::test_linear_jacobian(const P3D &pLocal, int qIdx) {
    Matrix dpdq_analytic, dpdq_estimated;
    compute_dpdq(pLocal, qIdx, dpdq_analytic);
    estimate_linear_jacobian(pLocal, qIdx, dpdq_estimated);

    bool error = false;

    for (int i = 0; i < dpdq_analytic.rows(); i++)
        for (int j = 0; j < dpdq_analytic.cols(); j++) {
            double err = dpdq_analytic(i, j) - dpdq_estimated(i, j);
            if (fabs(err) > 0.0001) {
                std::cout << "Error at: " << i << " " << j
                          << ": analytic: " << dpdq_analytic(i, j)
                          << " estimated " << dpdq_estimated(i, j)
                          << " error: " << err << std::endl;
                error = true;
            }
        }

    return !error;
}

bool GeneralizedCoordinatesRepresentation::test_linear_jacobian(const V3D &vLocal, int qIdx) {
    Matrix dvdq_analytic, dvdq_estimated;
    compute_dvdq(vLocal, qIdx, dvdq_analytic);
    estimate_linear_jacobian(vLocal, qIdx, dvdq_estimated);

    bool error = false;

    for (int i = 0; i < dvdq_analytic.rows(); i++)
        for (int j = 0; j < dvdq_analytic.cols(); j++) {
            double err = dvdq_analytic(i, j) - dvdq_estimated(i, j);
            if (fabs(err) > 0.0001) {
                std::cout << "Error at: " << i << " " << j
                          << ": analytic: " << dvdq_analytic(i, j)
                          << " estimated " << dvdq_estimated(i, j)
                          << " error: " << err << std::endl;
                error = true;
            }
        }

    return !error;
}

// returns the world-coordinates velocity for point pLocal, which is specified
// in the local coordinates of the child link of qIdx (relative to its COM)
V3D GeneralizedCoordinatesRepresentation::pDot(const P3D &pLocal, int qIdx) {
    // pDot(q) = dp/dq * qDot
    Matrix dpdq;
    Vector res;

    compute_dpdq(pLocal, qIdx, dpdq);
    res = dpdq * qDot;

    return V3D(res[0], res[1], res[2]);
}

// returns the world-coordinates angular velocity of the child link of qIdx
V3D GeneralizedCoordinatesRepresentation::RDot(int qIdx) {
    // w(q) = dR/dq * qDot
    Matrix dRdq;
    Vector res;

    compute_dRdq(qIdx, dRdq);
    res = dRdq * qDot;

    return V3D(res[0], res[1], res[2]);
}

// computes the matrix that tells you how the jacobian dp/dq changes with
// respect to q_i. Returns true if it contains non-zero elements, false
// otherwise
bool GeneralizedCoordinatesRepresentation::compute_ddpdq_dqi(
        const P3D &pLocal, int qIdx, Matrix &ddpdq_dqi, int q_i) {
    resize(ddpdq_dqi, 3, (int)q.size());

    int startIndex = qIdx;

    // if q_i is not one of the ancestors of rb, then it means the jacobian dpdq
    // does not depent on it, so check first...
    int qIndex = startIndex;
    bool isAncestor = false;
    while (qIndex > 2) {
        if (q_i == qIndex) {
            isAncestor = true;
            break;
        }
        qIndex = getParentQIdxOf(qIndex);
    }

    if (!isAncestor) return false;

    // input is valid, so we must compute this derivative...
    int loopIndex = startIndex;
    // 2 here is the index of the first translational DOF of the root
    while (loopIndex > 2) {
        int qIndex = startIndex;

        int stopIndex = loopIndex;
        if (q_i > loopIndex) stopIndex = q_i;

        //this is the point in the local coordinate frame of rb
        V3D vLocal = qOffsetToChildLinkCOM[qIdx] + V3D(P3D(), pLocal);

        while (qIndex > stopIndex) {
            //we need to compute the coordinates of pLocal in the frame of the parent DOF
            vLocal = qOffsetToParent[qIndex] + Quaternion(q[qIndex], getQAxis(qIndex)) * vLocal;

            //and repeat... compute the local coordinates of the point in the frame of the parent after joint rotation is applied
            qIndex = getParentQIdxOf(qIndex);
        }

        //make sure the offset vector we're looking at here is from the point around which the link rotates to the point of interest
        V3D offset = vLocal;

        while (qIndex > 2) {
            if (qIndex == loopIndex) offset = getQAxis(qIndex).cross(offset);
            if (qIndex == q_i) offset = getQAxis(qIndex).cross(offset);

            offset = Quaternion(q[qIndex], getQAxis(qIndex)) * offset;
            qIndex = getParentQIdxOf(qIndex);
        }

        for (int i = 0; i < 3; i++) ddpdq_dqi(i, loopIndex) = offset[i];

        loopIndex = getParentQIdxOf(loopIndex);
    }

    return true;
}

// computes the matrix that tells you how the jacobian dv/dq changes with
// respect to q_i. Returns true if it contains non-zero elements, false
// otherwise
bool GeneralizedCoordinatesRepresentation::compute_ddvdq_dqi(
        const V3D &vLocal, int qIdx, Matrix &ddvdq_dqi, int q_i) {
    resize(ddvdq_dqi, 3, (int)q.size());

    int startIndex = qIdx;

    // if q_i is not one of the ancestors of rb, then it means the jacobian dpdq
    // does not depent on it, so check first...
    int qIndex = startIndex;
    bool isAncestor = false;
    while (qIndex > 2) {
        if (q_i == qIndex) {
            isAncestor = true;
            break;
        }
        qIndex = getParentQIdxOf(qIndex);
    }

    if (!isAncestor) return false;

    // input is valid, so we must compute this derivative...
    int loopIndex = startIndex;
    // 2 here is the index of the first translational DOF of the root
    while (loopIndex > 2) {
        V3D theVector = vLocal;
        int qIndex = startIndex;

        int stopIndex = loopIndex;
        if (q_i > loopIndex) stopIndex = q_i;

        while (qIndex > stopIndex) {
            theVector = Quaternion(q[qIndex], getQAxis(qIndex)) * theVector;
            qIndex = getParentQIdxOf(qIndex);
        }

        while (qIndex > 2) {
            if (qIndex == loopIndex)
                theVector = getQAxis(qIndex).cross(theVector);
            if (qIndex == q_i) theVector = getQAxis(qIndex).cross(theVector);

            theVector = Quaternion(q[qIndex], getQAxis(qIndex)) * theVector;
            qIndex = getParentQIdxOf(qIndex);
        }

        for (int i = 0; i < 3; i++) ddvdq_dqi(i, loopIndex) = theVector[i];

        loopIndex = getParentQIdxOf(loopIndex);
    }

    return true;
}

// estimates the change of dp/dq with respect to q_i
void GeneralizedCoordinatesRepresentation::estimate_ddpdq_dqi(
        const P3D &pLocal, int qIdx, Matrix &ddpdq_dqi, int q_i) {
    resize(ddpdq_dqi, 3, (int)q.size());
    Matrix dpdq_p, dpdq_m;
    dpdq_p = ddpdq_dqi;
    dpdq_m = ddpdq_dqi;

    double val = q[q_i];
    double h = 0.0001;

    q[q_i] = val + h;
    compute_dpdq(pLocal, qIdx, dpdq_p);

    q[q_i] = val - h;
    compute_dpdq(pLocal, qIdx, dpdq_m);

    q[q_i] = val;

    ddpdq_dqi = (dpdq_p - dpdq_m) / (2 * h);
}

// estimates the change of dv/dq with respect to q_i
void GeneralizedCoordinatesRepresentation::estimate_ddvdq_dqi(
        const V3D &vLocal, int qIdx, Matrix &ddvdq_dqi, int q_i) {
    resize(ddvdq_dqi, 3, (int)q.size());
    Matrix dvdq_p, dvdq_m;
    dvdq_p = ddvdq_dqi;
    dvdq_m = ddvdq_dqi;

    double val = q[q_i];
    double h = 0.0001;

    q[q_i] = val + h;
    compute_dvdq(vLocal, qIdx, dvdq_p);

    q[q_i] = val - h;
    compute_dvdq(vLocal, qIdx, dvdq_m);

    q[q_i] = val;

    ddvdq_dqi = (dvdq_p - dvdq_m) / (2 * h);
}

bool GeneralizedCoordinatesRepresentation::
    test_linear_jacobian_derivatives(const P3D &pLocal, int qIdx) {
    Matrix dpdq_analytic, dpdq_estimated;

    bool error = false;

    for (int k = 0; k < (int)q.size(); k++) {
        compute_ddpdq_dqi(pLocal, qIdx, dpdq_analytic, k);
        estimate_ddpdq_dqi(pLocal, qIdx, dpdq_estimated, k);

        //		dpdq_analytic.printMatrix("out\\dpdq.mat");

        for (int i = 0; i < dpdq_analytic.rows(); i++)
            for (int j = 0; j < dpdq_analytic.cols(); j++) {
                double err = dpdq_analytic(i, j) - dpdq_estimated(i, j);
                if (fabs(err) > 0.0001) {
                    std::cout << "Error when computing ddpdq_dq " << k
                              << " at: " << i << " " << j
                              << ": analytic: " << dpdq_analytic(i, j)
                              << " estimated " << dpdq_estimated(i, j)
                              << " error: " << err << std::endl;
                    error = true;
                }
            }
    }

    return !error;
}

bool GeneralizedCoordinatesRepresentation::
    test_linear_jacobian_derivatives(const V3D &vLocal, int qIdx) {
    Matrix dvdq_analytic, dvdq_estimated;

    bool error = false;

    for (int k = 0; k < (int)q.size(); k++) {
        compute_ddvdq_dqi(vLocal, qIdx, dvdq_analytic, k);
        estimate_ddvdq_dqi(vLocal, qIdx, dvdq_estimated, k);

        //		dpdq_analytic.printMatrix("out\\dpdq.mat");

        for (int i = 0; i < dvdq_analytic.rows(); i++)
            for (int j = 0; j < dvdq_analytic.cols(); j++) {
                double err = dvdq_analytic(i, j) - dvdq_estimated(i, j);
                if (fabs(err) > 0.0001) {
                    std::cout << "Error when computing ddvdq_dq " << k
                              << " at: " << i << " " << j
                              << ": analytic: " << dvdq_analytic(i, j)
                              << " estimated " << dvdq_estimated(i, j)
                              << " error: " << err << std::endl;
                    error = true;
                }
            }
    }

    return !error;
}

// computes d(Jw)/dqi, where Jw is dR/dq
bool GeneralizedCoordinatesRepresentation::compute_ddRdq_dqi(int qIdx, Matrix &ddRdqdqi, int q_i) {
    resize(ddRdqdqi, 3, (int)q.size());
    if (q_i < 3) return false;

    int startIndex = qIdx;

    std::vector<int> qIndices;

    // if q_i is not one of the ancestors of rb, then it means the jacobian dpdq
    // does not depent on it, so check first...
    int qIndex = startIndex;
    bool isAncestor = false;

    Quaternion Q_qi = Quaternion::Identity();

    while (qIndex > 2) {
        qIndices.push_back(qIndex);

        if (q_i == qIndex) {
            isAncestor = true;
        }
        if (isAncestor) Q_qi = Quaternion(q[qIndex], getQAxis(qIndex)) * Q_qi;
        qIndex = getParentQIdxOf(qIndex);
    }

    if (!isAncestor) return false;

    V3D q_i_rotAxisWorld = Q_qi * getQAxis(q_i);

    Quaternion Q = Quaternion::Identity();
    // q rotates the rigid body (grand...child) about its world coordinate
    // axis... these will be the entries of dR/dq...
    while (qIndices.size() > 0) {
        int index = qIndices.back();
        Q = Q * Quaternion(q[qIndices.back()], getQAxis(qIndices.back()));

        if (index > q_i) {
            V3D derivativeVal = q_i_rotAxisWorld.cross(Q * getQAxis(index));

            ddRdqdqi(0, index) = derivativeVal(0);
            ddRdqdqi(1, index) = derivativeVal(1);
            ddRdqdqi(2, index) = derivativeVal(2);
        }
        qIndices.pop_back();
    }

    /*

            //input is valid, so we must compute this derivative...
            int loopIndex = startIndex;
            //2 here is the index of the first translational DOF of the root
            while (loopIndex > q_i) {
    //		Logger::consolePrint("NEW: adding contribution for index: %d\n",
    loopIndex); V3D derivativeVal =
    getWorldCoordsAxisForQ(q_i).cross(getWorldCoordsAxisForQ(loopIndex));
                    ddRdqdqi(0, loopIndex) = derivativeVal(0);
                    ddRdqdqi(1, loopIndex) = derivativeVal(1);
                    ddRdqdqi(2, loopIndex) = derivativeVal(2);

                    loopIndex = qParentIndex[loopIndex];
            }
    */
    return true;
}

// estimates the change of angular jacobian with respect to q_i using finite
// differences
void GeneralizedCoordinatesRepresentation::estimate_ddRdq_dqi(
        int qIdx, Matrix &ddRdq_dqi, int q_i) {
    resize(ddRdq_dqi, 3, (int)q.size());
    Matrix dRdq_p, dRdq_m;
    dRdq_p = ddRdq_dqi;
    dRdq_m = ddRdq_dqi;

    double val = q[q_i];
    double h = 0.0001;

    q[q_i] = val + h;
    compute_dRdq(qIdx, dRdq_p);

    q[q_i] = val - h;
    compute_dRdq(qIdx, dRdq_m);

    q[q_i] = val;

    ddRdq_dqi = (dRdq_p - dRdq_m) / (2 * h);
}

bool GeneralizedCoordinatesRepresentation::test_angular_jacobian_derivatives(int qIdx) {
    Matrix ddRdqdq_analytic, ddRdqdq_estimated;

    bool error = false;

    for (int k = 0; k < (int)q.size(); k++) {
        compute_ddRdq_dqi(qIdx, ddRdqdq_analytic, k);
        estimate_ddRdq_dqi(qIdx, ddRdqdq_estimated, k);

        // print("../out/angular_jacobian_dqi_analytic.mat",
        // ddRdqdq_analytic);
        // print("../out/angular_jacobian_dqi_estimated.mat",
        // ddRdqdq_estimated);

        for (int i = 0; i < ddRdqdq_analytic.rows(); i++)
            for (int j = 0; j < ddRdqdq_analytic.cols(); j++) {
                double err = ddRdqdq_analytic(i, j) - ddRdqdq_estimated(i, j);
                if (fabs(err) > 0.0001) {
                    std::cout << "Error when computing ddRdqdq " << k
                              << " at: " << i << " " << j
                              << ": analytic: " << ddRdqdq_analytic(i, j)
                              << " estimated " << ddRdqdq_estimated(i, j)
                              << " error: " << err << std::endl;
                    error = true;
                    //					exit(0);
                }
            }
    }

    return !error;
}

void testGeneralizedCoordinateRepresentation(ArticulatedFigure *af) {
    std::cout << "testing robot generalized coordinates representation..."
              << std::endl;

    // make sure we project errors introduced by physics engine (i.e. hinge
    // joints not rotating only about their axis)
    af->fixArticulationConstraints();

    // test out projections between robot state and generalized coordinates...
    AFState robotState1(af);

    GeneralizedCoordinatesRepresentation gcrrNew(af);

    AFState robotState2 = gcrrNew.getAFState();

    if (!(robotState1 == robotState2)) {
        std::cout << "TESTING GENERALIZED COORDINATES: robot state is not the "
                     "same after projection..."
                  << std::endl;
    }

    // test forward kinematics (world pos of points on RBs, velocity of points,
    // etc)...
    for (int i = 0; i < af->getJointCount(); i++) {
//        std::cout << "Testing joint " << af->joints[i]->name << std::endl;

        P3D point = P3D(0, 0, 0) + V3D(getRandomNumberIn01Range(), getRandomNumberIn01Range(), getRandomNumberIn01Range()) * 0.2;
        P3D wc1 = af->joints[i]->child->getWorldCoordinates(point);
        P3D wc2 = gcrrNew.p_of_q(point, 6 + af->idx(af->getParentOf(af->joints[i]->child)));
        if (V3D(wc1, wc2).norm() > 1e-8)
            std::cout << "TESTING GENERALIZED COORDINATES: world coordinates "
                         "of point on rigid body" << af->joints[i]->child->name << "do not match up... error: "
                      << V3D(wc1, wc2).norm() << std::endl;

        wc1 = af->joints[i]->parent->getWorldCoordinates(point);
        wc2 = gcrrNew.p_of_q(point, 6 + af->idx(af->getParentOf(af->joints[i]->parent)));
        if (V3D(wc1, wc2).norm() > 1e-8)
            std::cout << "TESTING GENERALIZED COORDINATES: world coordinates "
                         "of point on rigid body " << af->joints[i]->parent->name << "do not match up... error: "
                      << V3D(wc1, wc2).norm() << std::endl;

        V3D vec = V3D(getRandomNumberIn01Range(), getRandomNumberIn01Range(), getRandomNumberIn01Range()) * 0.2;
        V3D vc1 = af->joints[i]->child->getWorldCoordinates(vec);
        V3D vc2 = gcrrNew.v_of_q(vec, 6 + af->idx(af->getParentOf(af->joints[i]->child)));

        if ((vc1 - vc2).norm() > 1e-8)
            std::cout << "TESTING GENERALIZED COORDINATES: world coordinates "
                         "of vectors on rigid body " << af->joints[i]->child->name << "do not match up... error: "
                      << (vc1 - vc2).norm() << std::endl;

        vc1 = af->joints[i]->parent->getWorldCoordinates(vec);
        vc2 = gcrrNew.v_of_q(vec, 6 + af->idx(af->getParentOf(af->joints[i]->parent)));

        if ((vc1 - vc2).norm() > 1e-8)
            std::cout << "TESTING GENERALIZED COORDINATES: world coordinates "
                         "of vectors on rigid body " << af->joints[i]->parent->name << " do not match up... error: "
                      << (vc1 - vc2).norm() << std::endl;

        V3D wv1 = af->joints[i]->child->getVelocityForLocalCoordsPoint(point);
        V3D wv2 = gcrrNew.pDot(point, 6 + af->idx(af->getParentOf(af->joints[i]->child)));
        if ((wv1 - wv2).norm() > 1e-8)
            std::cout << "TESTING GENERALIZED COORDINATES: velocities of point "
                         "on rigid body do not match up... error: "
                      << (wv1 - wv2).norm() << std::endl;

        wv1 = af->joints[i]->parent->velocity;
        wv2 = gcrrNew.pDot(P3D(0, 0, 0), 6 + af->idx(af->getParentOf(af->joints[i]->parent)));
        if ((wv1 - wv2).norm() > 1e-8)
            std::cout << "TESTING GENERALIZED COORDINATES: velocities of rigid "
                         "body do not match up... error: "
                      << (wv1 - wv2).norm() << std::endl;

        wv1 = af->joints[i]->child->angularVelocity;
        wv2 = gcrrNew.RDot(6 + af->idx(af->getParentOf(af->joints[i]->child)));
        if ((wv1 - wv2).norm() > 1e-8)
            std::cout << "TESTING GENERALIZED COORDINATES: angular velocities "
                         "of rigid body do not match up... error: "
                      << (wv1 - wv2).norm() << std::endl;

        Quaternion r1 = af->joints[i]->child->orientation;
        Quaternion r2 = gcrrNew.R_of_q(6 + af->idx(af->getParentOf(af->joints[i]->child)));
        if ((r1.inverse() * r2).vec().norm() > 1e-8) {
            std::cout << "TESTING GENERALIZED COORDINATES: orientations of "
                         "rigid body do not match up... error: "
                      << r1.w() << " / " << r2.w() << std::endl;
        }

        if (!gcrrNew.test_linear_jacobian(point, 6 + af->idx(af->getParentOf(af->joints[i]->child))))
            std::cout << "TESTING GENERALIZED COORDINATES: linear jacobian "
                         "(P3D) does not match FD..."
                      << std::endl;

        if (!gcrrNew.test_linear_jacobian(point, 6 + af->idx(af->getParentOf(af->joints[i]->parent))))
            std::cout << "TESTING GENERALIZED COORDINATES: linear jacobian "
                         "(P3D) does not match FD..."
                      << std::endl;

        if (!gcrrNew.test_linear_jacobian(vec, 6 + af->idx(af->getParentOf(af->joints[i]->child))))
            std::cout << "TESTING GENERALIZED COORDINATES: linear jacobian "
                         "(V3D) does not match FD..."
                      << std::endl;

        if (!gcrrNew.test_linear_jacobian(vec, 6 + af->idx(af->getParentOf(af->joints[i]->parent))))
            std::cout << "TESTING GENERALIZED COORDINATES: linear jacobian "
                         "(V3D) does not match FD..."
                      << std::endl;

        if (!gcrrNew.test_angular_jacobian(6 + af->idx(af->getParentOf(af->joints[i]->child))))
            std::cout << "TESTING GENERALIZED COORDINATES: angular jacobian "
                         "does not match FD..."
                      << std::endl;

        if (!gcrrNew.test_linear_jacobian_derivatives(
                point, 6 + af->idx(af->getParentOf(af->joints[i]->child))))
            std::cout << "TESTING GENERALIZED COORDINATES: linear jacobian "
                         "derivatives do not match FD..."
                      << std::endl;

        if (!gcrrNew.test_linear_jacobian_derivatives(
                vec, 6 + af->idx(af->getParentOf(af->joints[i]->child))))
            std::cout << "TESTING GENERALIZED COORDINATES: linear jacobian "
                         "derivatives do not match FD..."
                      << std::endl;

        if (!gcrrNew.test_linear_jacobian_derivatives(
                point, 6 + af->idx(af->getParentOf(af->joints[i]->parent))))
            std::cout << "TESTING GENERALIZED COORDINATES: linear jacobian "
                         "derivatives do not match FD..."
                      << std::endl;

        if (!gcrrNew.test_linear_jacobian_derivatives(
                vec, 6 + af->idx(af->getParentOf(af->joints[i]->parent))))
            std::cout << "TESTING GENERALIZED COORDINATES: linear jacobian "
                         "derivatives do not match FD..."
                      << std::endl;

        if (!gcrrNew.test_angular_jacobian_derivatives(
                6 + af->idx(af->getParentOf(af->joints[i]->child))))
            std::cout << "TESTING GENERALIZED COORDINATES: angular jacobian "
                         "derivatives do not match FD..."
                      << std::endl;

    }
}

