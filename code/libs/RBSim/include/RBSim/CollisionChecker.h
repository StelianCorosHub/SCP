#pragma once

#include <utils/geoms.h>
#include <utils/mathUtils.h>
#include <utils/utils.h>
#include <imgui.h>

class ContactPoint {
public:
    P3D p;
    V3D n = V3D(0,1,0);
    double d = 0;

    ContactPoint() {}

    ContactPoint(const P3D& p, const V3D& n, double d) {
        if (d > 0)
            Logger::consolePrint("Warning: adding contact point with positive penetration depth...");

        this->p = p;
        this->n = n;
        this->d = d;
    }

    ContactPoint(const P3D& p) {
        this->p = p;
    }
};

class DebugContactPoint {
public:
    P3D p;
    V3D color;

    DebugContactPoint(const P3D& p, const V3D& color = V3D(1,0,0)) {
        this->p = p;
        this->color = color;
    }

};

// find all the intersection points between the 2D rectangle with vertices
// at (+/-h[0],+/-h[1]) and the 2D quadrilateral with vertices (p[0],p[1]),
// (p[2],p[3]),(p[4],p[5]),(p[6],p[7]).
//
// the intersection points are returned as x,y pairs in the 'ret' array.
// the number of intersection points is returned by the function (this will
// be in the range 0 to 8).

inline int intersectRectQuad2 (double h[2], double p[8], double ret[16]) {
// q (and r) contain nq (and nr) coordinate points for the current (and chopped) polygons
    int nq=4, nr=0;
    double buffer[16];
    double *q = p;
    double *r = ret;

    // direction notation: xy[0] = x axis, xy[1] = y axis
    for (int dir=0; dir <= 1; dir++) {
        for (int sign=-1; sign <= 1; sign += 2) {
            // chop q along the line xy[dir] = sign*h[dir]
            double *pq = q;
            double *pr = r;
            nr = 0;
            for (int i=nq; i > 0; i--) {
                // go through all points in q and all lines between adjacent points
                if (sign*pq[dir] < h[dir]) {
                    // this point is inside the chopping line
                    pr[0] = pq[0];
                    pr[1] = pq[1];
                    pr += 2;
                    nr++;
                    if (nr & 8) {
                        q = r;
                        goto done;
                    }
                }
                double *nextq = (i > 1) ? pq+2 : q;
                if ((sign*pq[dir] < h[dir]) ^ (sign*nextq[dir] < h[dir])) {
                    // this line crosses the chopping line
                    pr[1-dir] = pq[1-dir] + (nextq[1-dir]-pq[1-dir]) / (nextq[dir]-pq[dir]) * (sign*h[dir]-pq[dir]);
                    pr[dir] = sign*h[dir];
                    pr += 2;
                    nr++;
                    if (nr & 8) {
                        q = r;
                        goto done;
                    }
                }
                pq += 2;
            }
            q = r;
            r = (q==ret) ? buffer : ret;
            nq = nr;
        }
    }
done:
    if (q != ret) memcpy (ret,q,nr*2*sizeof(double));
        return nr;
}

// Given two boxes (pA,qA,dA) and (pB,qB,dB), where p is the position of
// the center of a box, q represents its orientation, and d is a vector
// that stores its dimension (length, width, height), generate a list
// of contact points between them. For each contact point, we will return the
// normal (pointing away from box A and towards box b) and the pentration depth.

inline Array<ContactPoint> collideBoxWithBox ( const P3D& pA, const Quaternion& qA, V3D dA,
                                const P3D& pB, const Quaternion& qB, V3D dB,
                                Array<DebugContactPoint> *dCPs = nullptr) {

    Array<ContactPoint> CPs;

// we need to find a separating plane between the two boxes. If one is found, then
// there is no intersection. Otherwise, we'll need to go ahead and compute intersection
// points.

// There are 15 separating planes to check, and the checks can be done quite efficiently
// based only on their normal vector:
// https://jkh.me/files/tutorials/Separating%20Axis%20Theorem%20for%20Oriented%20Bounding%20Boxes.pdf

    dA /= 2;
    dB /= 2;

    //vector from centers of box 1 to box 2
    V3D T(pA, pB);
    //in these rotation matrices, the columns represent the x-, y- and z- axes of the boxes in world coordinates
    Matrix3x3 RA(qA), RB(qB);
    double rA = MAX(dA.x(), dA.y()); rA = MAX(rA, dA.z());
    double rB = MAX(dB.x(), dB.y()); rB = MAX(rB, dB.z());

    if (T.norm() > dA.norm() + dB.norm())
        return CPs;

#define TEST_SEPARATING_PLANE(N)                                                    \
    {                                                                               \
        cIdxTmp++;                                                                  \
        V3D n = V3D(N);                                                             \
        double nLen = n.norm();                                                     \
        if (nLen > 1e-5) {                                                          \
            n /= nLen;                                                              \
            double rhs = fabs(T.dot(n));                                            \
            double lhs = 0;                                                         \
            for (int i = 0; i < 3; i++) {                                           \
                lhs += fabs(dA(i) * RA.col(i).dot(n));                              \
                lhs += fabs(dB(i) * RB.col(i).dot(n));                              \
            }                                                                       \
            double pDepth = rhs - lhs;                                              \
            if (pDepth > 0)                                                         \
                return CPs;                                                         \
            if (pDepth * safetyFactor > pVal) {                                     \
                separationPlaneN = n;                                               \
                if (T.dot(n) < 0) {                                                 \
                    separationPlaneN *= -1;                                         \
                }                                                                   \
                pVal = pDepth;                                                      \
                codeIdx = cIdxTmp;                                                  \
            }                                                                       \
        }                                                                           \
    }

    int cIdxTmp = -1;

    //we will want to keep track of the plane that is closest to separating the two boxes
    double pVal = -DBL_MAX;
    //normal for the plane that separates the two boxes the most - pointing mostly from box a towards box b
    V3D separationPlaneN;
    int codeIdx = 0;

    double safetyFactor = 1.0;

    //here are the 15 plane normals that can define a separation plane...
    TEST_SEPARATING_PLANE(RA.col(0));
    TEST_SEPARATING_PLANE(RA.col(1));
    TEST_SEPARATING_PLANE(RA.col(2));

    TEST_SEPARATING_PLANE(RB.col(0));
    TEST_SEPARATING_PLANE(RB.col(1));
    TEST_SEPARATING_PLANE(RB.col(2));

    //edge-edge intersection tests are more finicky... if we have a good-enough face-something intersection test
    //we should be prioritizing that one...
    safetyFactor = 5.0;

    TEST_SEPARATING_PLANE(RA.col(0).cross(RB.col(0)));
    TEST_SEPARATING_PLANE(RA.col(0).cross(RB.col(1)));
    TEST_SEPARATING_PLANE(RA.col(0).cross(RB.col(2)));

    TEST_SEPARATING_PLANE(RA.col(1).cross(RB.col(0)));
    TEST_SEPARATING_PLANE(RA.col(1).cross(RB.col(1)));
    TEST_SEPARATING_PLANE(RA.col(1).cross(RB.col(2)));

    TEST_SEPARATING_PLANE(RA.col(2).cross(RB.col(0)));
    TEST_SEPARATING_PLANE(RA.col(2).cross(RB.col(1)));
    TEST_SEPARATING_PLANE(RA.col(2).cross(RB.col(2)));

//    Logger::consolePrint("collision code: %d\n", codeIdx);

    //depending on which separation plane we've got we'll be computing collision points in different ways...

    if (codeIdx >= 6) {
        // an edge from box 1 touches an edge from box 2.

        // the code cIdx and normal will tell us exactly which edges these are...
        int bAEdgeIdx = (codeIdx - 6) / 3;
        int bBEdgeIdx = (codeIdx - 6) % 3;

        V3D bAEdgeV = V3D(RA.col(bAEdgeIdx));
        V3D bBEdgeV = V3D(RB.col(bBEdgeIdx));

        //now, choose a point on the correct edge for each of the boxes
        P3D pointOnbA = pA;
        for (int i = 0; i < 3; i++)
            if (i != bAEdgeIdx)
                pointOnbA += RA.col(i) * dA[i] * ((separationPlaneN.dot(RA.col(i)) < 0) ? -1 : 1);

        P3D pointOnbB = pB;
        for (int i = 0; i < 3; i++)
            if (i != bBEdgeIdx)
                pointOnbB += RB.col(i) * dB[i] * ((separationPlaneN.dot(RB.col(i)) > 0) ? -1 : 1);

        Segment s = Line(pointOnbA, bAEdgeV).getClosestSegmentToLine(Line(pointOnbB, bBEdgeV));

        if (dCPs) {
            dCPs->push_back(DebugContactPoint(pointOnbA + RA.col(bAEdgeIdx) * dA[bAEdgeIdx] * 1, V3D(1,0,1)));
            dCPs->push_back(DebugContactPoint(pointOnbA + RA.col(bAEdgeIdx) * dA[bAEdgeIdx] * -1, V3D(1,0,1)));
            dCPs->push_back(DebugContactPoint(pointOnbB + RB.col(bBEdgeIdx) * dB[bBEdgeIdx] * 1, V3D(1,1,0)));
            dCPs->push_back(DebugContactPoint(pointOnbB + RB.col(bBEdgeIdx) * dB[bBEdgeIdx] * -1, V3D(1,1,0)));
        }

        CPs.push_back(ContactPoint(s.b, separationPlaneN, pVal));
//        CPs.push_back(ContactPoint(s.a, nSel, pVal));

        return CPs;
    }

    // if we got here, then we have a face-something intersection (because the separating
    // plane is one of the faces of either box a or box b. define face 'a' to be the reference
    // face (i.e. the normal vector is perpendicular to this) and face 'b' to be
    // the incident face (the closest face of the other box). For convenience, R_ref, p_ref and S_ref
    // are the orientation / position / dimension of the box the reference face lies on, and ditto
    // for R_inc, p_inc and S_inc.

    Matrix3x3 *R_ref, *R_inc;
    const P3D *p_ref, *p_inc;
    const V3D *S_ref, *S_inc;

    R_ref = (codeIdx < 3 ) ? &RA : &RB;
    R_inc = (codeIdx < 3 ) ? &RB : &RA;
    p_ref = (codeIdx < 3 ) ? &pA : &pB;
    p_inc = (codeIdx < 3 ) ? &pB : &pA;
    S_ref = (codeIdx < 3 ) ? &dA : &dB;
    S_inc = (codeIdx < 3 ) ? &dB : &dA;
    if (codeIdx >=3) separationPlaneN *= -1;

    //these are the indices of the normal of the reference face, as well as the two axes spanning the reference face
    int refFaceNIdx = codeIdx % 3;
    int refFaceAxis1Idx = (refFaceNIdx == 0) ? 1 : 0;
    int refFaceAxis2Idx = refFaceAxis1Idx + 1;
    if (refFaceNIdx == refFaceAxis2Idx) refFaceAxis2Idx++;

    //these are in local coordinates
    V3D refFaceAxis1, refFaceAxis2;
    refFaceAxis1[refFaceAxis1Idx] = 1.0;
    refFaceAxis2[refFaceAxis2Idx] = 1.0;

    //this is the center of the reference face, in world coordinates
    P3D centerRefFace = *p_ref + separationPlaneN * (*S_ref)[refFaceNIdx];

    if (dCPs) {
        dCPs->push_back(DebugContactPoint(centerRefFace, V3D(1, 0, 1)));
        dCPs->push_back(DebugContactPoint(centerRefFace + (*R_ref) * (refFaceAxis1 * (*S_ref)[refFaceAxis1Idx] + refFaceAxis2 * (*S_ref)[refFaceAxis2Idx]), V3D(1, 0.5, 1)));
        dCPs->push_back(DebugContactPoint(centerRefFace + (*R_ref) * (refFaceAxis1 * (*S_ref)[refFaceAxis1Idx] - refFaceAxis2 * (*S_ref)[refFaceAxis2Idx]), V3D(1, 0.5, 1)));
        dCPs->push_back(DebugContactPoint(centerRefFace + (*R_ref) * (-refFaceAxis1 * (*S_ref)[refFaceAxis1Idx] + refFaceAxis2 * (*S_ref)[refFaceAxis2Idx]), V3D(1, 0.5, 1)));
        dCPs->push_back(DebugContactPoint(centerRefFace + (*R_ref) * (-refFaceAxis1 * (*S_ref)[refFaceAxis1Idx] - refFaceAxis2 * (*S_ref)[refFaceAxis2Idx]), V3D(1, 0.5, 1)));
    }

    //intersections will be computed in the plane of this reference face, where dimensions are
    double refFaceRect[2];
    refFaceRect[0] = (*S_ref)[refFaceAxis1Idx];
    refFaceRect[1] = (*S_ref)[refFaceAxis2Idx];

    //now we need to compute the incident face. We will get the one whose normal is most aligned with the normal of the reference face
    double dotN0 = fabs(separationPlaneN.dot(R_inc->col(0)));
    double dotN1 = fabs(separationPlaneN.dot(R_inc->col(1)));
    double dotN2 = fabs(separationPlaneN.dot(R_inc->col(2)));

    int incFaceNIdx = 0;
    if (dotN1 >= dotN0 && dotN1 >= dotN2) incFaceNIdx = 1;
    if (dotN2 >= dotN0 && dotN2 >= dotN1) incFaceNIdx = 2;

    int incFaceAxis1Idx = (incFaceNIdx == 0) ? 1 : 0;
    int incFaceAxis2Idx = incFaceAxis1Idx + 1;
    if (incFaceNIdx == incFaceAxis2Idx) incFaceAxis2Idx++;

    V3D incFaceN(R_inc->col(incFaceNIdx) * ((separationPlaneN.dot(R_inc->col(incFaceNIdx)) > 0) ? -1 : 1));

    P3D centerIncFace = *p_inc + incFaceN * (*S_inc)[incFaceNIdx];

    if (dCPs) {
        dCPs->push_back(DebugContactPoint(centerIncFace, V3D(1, 1, 0)));

        dCPs->push_back(DebugContactPoint(centerIncFace + R_inc->col(incFaceAxis1Idx) * (*S_inc)[incFaceAxis1Idx] + R_inc->col(incFaceAxis2Idx) * (*S_inc)[incFaceAxis2Idx], V3D(1, 1, 0.5)));
        dCPs->push_back(DebugContactPoint(centerIncFace + R_inc->col(incFaceAxis1Idx) * (*S_inc)[incFaceAxis1Idx] - R_inc->col(incFaceAxis2Idx) * (*S_inc)[incFaceAxis2Idx], V3D(1, 1, 0.5)));
        dCPs->push_back(DebugContactPoint(centerIncFace - R_inc->col(incFaceAxis1Idx) * (*S_inc)[incFaceAxis1Idx] + R_inc->col(incFaceAxis2Idx) * (*S_inc)[incFaceAxis2Idx], V3D(1, 1, 0.5)));
        dCPs->push_back(DebugContactPoint(centerIncFace - R_inc->col(incFaceAxis1Idx) * (*S_inc)[incFaceAxis1Idx] - R_inc->col(incFaceAxis2Idx) * (*S_inc)[incFaceAxis2Idx], V3D(1, 1, 0.5)));
    }

    //we now need to project each of the corner points of the incident face onto the reference face (in a coordinate frame
    //centered at the centerRefFace with axes given by refFaceAxis1 and refFaceAxis2
    double incFaceProjected[8];     // 2D coordinates of corners of incident face projected onto reference face, (x,y) pairs
    int idx = 0;
    V3D offset;

#define ADD_INCIDENT_FACE_CORNER_TO_QUAD(x, y) {                                                                                                                                    \
        offset = V3D(centerRefFace, centerIncFace + x * R_inc->col(incFaceAxis1Idx) * (*S_inc)[incFaceAxis1Idx] + y * R_inc->col(incFaceAxis2Idx) * (*S_inc)[incFaceAxis2Idx]);     \
        incFaceProjected[idx++] = offset.dot(R_ref->col(refFaceAxis1Idx));                                                                                                          \
        incFaceProjected[idx++] = offset.dot(R_ref->col(refFaceAxis2Idx));                                                                                                          \
    }

//    this projection step here is broken... going from points in 3d to points in 2d to points in 3d

    ADD_INCIDENT_FACE_CORNER_TO_QUAD(-1, -1);
    ADD_INCIDENT_FACE_CORNER_TO_QUAD(-1, 1);
    ADD_INCIDENT_FACE_CORNER_TO_QUAD(1, 1);
    ADD_INCIDENT_FACE_CORNER_TO_QUAD(1, -1);

    //all right, now we compute the intersection between the two (planar) faces...
    // intersect the incident and reference faces
    double intersectionPts[16];
    int nres = intersectRectQuad2 (refFaceRect, incFaceProjected, intersectionPts);
    if (nres < 1) {
//        Logger::consolePrint("On no, we got here... 11\n");
        return CPs;          // this should not happen
    }

    //for each intersection point, project it back into 3d space, and check if it lies "below"
    // the plane of the reference face. If it does, it must mean it's an intersection point...
    for (int i = 0; i < nres; i++) {
        P3D p = centerRefFace + (*R_ref) * (refFaceAxis1 * intersectionPts[2 * i + 0] + refFaceAxis2 * intersectionPts[2 * i + 1]);

        if (dCPs)
            dCPs->push_back(DebugContactPoint(p, V3D(1, 1, 1)));

        //now, move this point up onto the plane of the incident face along the normal of the reference face
        double t = -incFaceN.dot(V3D(centerIncFace, p)) / incFaceN.dot(separationPlaneN);
        if (t < 0) {
            p += t * separationPlaneN;
            if (codeIdx >= 3)
                CPs.push_back(ContactPoint(p, -separationPlaneN, t));
            else
                CPs.push_back(ContactPoint(p, separationPlaneN, t));
        }
    }

    return CPs;
}

// Given a box (pB,qB,dB) and a sphere (pS, r), where p is the position of
// the center of the box/sphere, q is the orientation of the box, d is a vector
// that stores the dimension (length, width, height) of the box, and r is the
// radius of the sphere, generate a list of contact points between them. If a contact
// point is generated, the normal points away from the box towards the sphere and the
// pentration depth, a negative scalar, tells us how deep the collision is.
inline Array<ContactPoint> collideBoxWithSphere ( const P3D& pB, const Quaternion& qB, V3D dB,
                                               const P3D& pS, double r) {

    Array<ContactPoint> CPs;

    dB *= 0.5;

    //compute the position of the sphere in the coordinate frame of the box
    P3D pLocal = P3D() + qB.inverse() * V3D(pB, pS);

    if (pLocal.x >= -dB.x() && pLocal.x <= dB.x() &&
        pLocal.y >= -dB.y() && pLocal.y <= dB.y() &&
        pLocal.z >= -dB.z() && pLocal.z <= dB.z()) {

        //center of sphere is inside the box... now have to check each of the six planes of the box to see which one
        //the center point is closest to
        double dVal = -DBL_MAX;
        V3D n;

#define TEST_DIM(i, m){                                         \
                        double d = m * (pLocal[i] - dB[i] * m); \
                        if (d > dVal) {                         \
                            dVal = d;                           \
                            n = V3D(); n[i] = m;                \
                        }                                       \
                      }

        TEST_DIM(0, 1);
        TEST_DIM(0, -1);

        TEST_DIM(1, 1);
        TEST_DIM(1, -1);

        TEST_DIM(2, 1);
        TEST_DIM(2, -1);

        CPs.push_back(ContactPoint(pB + qB * V3D(P3D(), pLocal), qB * n, dVal));

    } else {
        //must check distance from the sphere to each of the planes of the box...

        P3D pProj = pLocal;
        for (int i = 0; i < 3; i++)
            clamp(pProj[i], -dB[i], dB[i]);

        V3D offset = V3D(pProj, pLocal);
        double d = offset.norm();

        if (d < r)
            CPs.push_back(ContactPoint(pB + qB * V3D(P3D(), pProj), qB * offset.unit(), d - r));
    }

    return CPs;
}

// Similar to collideBoxWithSphere, but normals point away from sphere, towards the box
inline Array<ContactPoint> collideSphereWithBox(const P3D& pS, double r,
                                                const P3D& pB, const Quaternion& qB, V3D dB) {

    Array<ContactPoint> CPs = collideBoxWithSphere(pB, qB, dB, pS, r);

    for (int i = 0; i < CPs.size(); i++)
        CPs[i].n *= -1;

    return CPs;
}

// collides the two spheres. The normal of the contacts points away from sphere A towards sphere B
inline Array<ContactPoint> collideSphereWithSphere(const P3D& pA, double rA, const P3D& pB, double rB) {
    V3D v(pA, pB);
    double d = v.norm() - rA - rB;
    v.normalize();

    Array<ContactPoint> CPs;

    if (d < 0)
        CPs.push_back(ContactPoint(pA + v * (rA + d / 2.0), v, d));

    return CPs;
}

// Collides box (pB,qB,dB) with the plane (pP,nP) and generates a list
// of contact points between them. For each contact point, we will return the
// normal (pointing away from the box) and the pentration depth.

inline Array<ContactPoint> collideBoxWithPlane (const P3D& pB, const Quaternion& qB, V3D dB,
                                               const P3D& pP, const V3D& nP) {

    Array<ContactPoint> CPs;

    dB /= 2;

    Plane p(pP, nP);

#define TEST_COL(pt) {                                          \
              P3D P = pt;                                       \
              double d = p.getSignedDistanceToPoint(P);         \
              if (d < 0)                                        \
                    CPs.push_back(ContactPoint(P, -nP, d));     \
            }

    //check each corner of the box and add it as a collision point if it lies below the plane...
    TEST_COL(pB + qB * (V3D(1,0,0) * dB[0] + V3D(0,1,0) * dB[1] + V3D(0,0,1) * dB[2]));
    TEST_COL(pB + qB * (V3D(1,0,0) * dB[0] + V3D(0,1,0) * dB[1] - V3D(0,0,1) * dB[2]));
    TEST_COL(pB + qB * (V3D(1,0,0) * dB[0] - V3D(0,1,0) * dB[1] + V3D(0,0,1) * dB[2]));
    TEST_COL(pB + qB * (V3D(1,0,0) * dB[0] - V3D(0,1,0) * dB[1] - V3D(0,0,1) * dB[2]));
    TEST_COL(pB - qB * (V3D(1,0,0) * dB[0] + V3D(0,1,0) * dB[1] + V3D(0,0,1) * dB[2]));
    TEST_COL(pB - qB * (V3D(1,0,0) * dB[0] + V3D(0,1,0) * dB[1] - V3D(0,0,1) * dB[2]));
    TEST_COL(pB - qB * (V3D(1,0,0) * dB[0] - V3D(0,1,0) * dB[1] + V3D(0,0,1) * dB[2]));
    TEST_COL(pB - qB * (V3D(1,0,0) * dB[0] - V3D(0,1,0) * dB[1] - V3D(0,0,1) * dB[2]));

    return CPs;
}

// Same as collideBoxWithPlane, but normals point towards the box
inline Array<ContactPoint> collidePlaneWithBox(const P3D& pP, const V3D& nP,
                                               const P3D& pB, const Quaternion& qB, V3D dB) {

    Array<ContactPoint> CPs = collideBoxWithPlane(pB, qB, dB, pP, nP);

    for (int i = 0; i < CPs.size(); i++)
        CPs[i].n *= -1;

    return CPs;
}

// Collides sphere (pS,r) with the plane (pP,nP) and generates a list
// of contact points between them. For each contact point, we will return the
// normal (pointing away from the sphere) and the pentration depth.
inline Array<ContactPoint> collideSphereWithPlane(const P3D& pS, double r,
                                                  const P3D& pP, const V3D& nP) {

    Array<ContactPoint> CPs;
    Plane p(pP, nP);

    double d = p.getSignedDistanceToPoint(pS);

    if (d < r)
        CPs.push_back(ContactPoint(pS - nP * r, -nP, d - r));

    return CPs;
}

// Same as collideSphereWithPlane, but normals point towards the sphere
inline Array<ContactPoint> collidePlaneWithSphere(const P3D& pP, const V3D& nP, const P3D& pS, double r) {

    Array<ContactPoint> CPs = collideSphereWithPlane(pS, r, pP, nP);

    for (int i = 0; i < CPs.size(); i++)
        CPs[i].n *= -1;

    return CPs;
}

inline Array<ContactPoint> collideSphereWithSphere(const Sphere& a, const Sphere& b) {
    return collideSphereWithSphere(a.p, a.r, b.p, b.r);
}

inline Array<ContactPoint> collideSphereWithBox(const Sphere& a, const Box3D& b) {
    return collideSphereWithBox(a.p, a.r, b.p, b.q, b.dims);
}

inline Array<ContactPoint> collideSphereWithPlane(const Sphere& a, const Plane& b) {
    return collideSphereWithPlane(a.p, a.r, b.p, b.n);
}

inline Array<ContactPoint> collideBoxWithSphere(const Box3D& a, const Sphere& b) {
    return collideBoxWithSphere(a.p, a.q, a.dims, b.p, b.r);
}

inline Array<ContactPoint> collideBoxWithBox(const Box3D& a, const Box3D& b) {
    return collideBoxWithBox(a.p, a.q, a.dims, b.p, b.q, b.dims);
}

inline Array<ContactPoint> collideBoxWithPlane(const Box3D& a, const Plane& b) {
    return collideBoxWithPlane(a.p, a.q, a.dims, b.p, b.n);
}

inline Array<ContactPoint> collidePlaneWithSphere(const Plane& a, const Sphere& b) {
    return collidePlaneWithSphere(a.p, a.n, b.p, b.r);
}

inline Array<ContactPoint> collidePlaneWithBox(const Plane& a, const Box3D& b) {
    return collidePlaneWithBox(a.p, a.n, b.p, b.q, b.dims);
}
