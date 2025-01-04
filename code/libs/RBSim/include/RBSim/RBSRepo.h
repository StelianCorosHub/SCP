#pragma once

#include <RBSim/RBJoint.h>
#include <RBSim/RB.h>
#include <RBSim/RBUtils.h>
#include <utils/mathUtils.h>
#include <utils/utils.h>

#include <RBSim/RBSRenderer.h>

#include <unordered_map>

/**
 * This class stores a list of rigid bodies and joints using smart pointers
 * and implements useful getters. It will be extended by robots, or physics
 * engines, for which collections of rigid bodies and joints are the foundation.
 */
class RBSRepo {
protected:
    // quick access to the indices of all RBs and joints in the
    // global lists managed by the articulated figure data structure
    std::unordered_map<pRB, int> rbIndices;
    std::unordered_map<pRBJ, int> jIndices;

public:

    RBSRepo() {}

    ~RBSRepo() {}

    Array<pRB> RBs;
    Array<pRBJ> joints;

    DrawingFlagContainer df;

    /**
     * This method returns the shared pointer reference to the rigid body with
     * the given name, or nullptr if it is not found
     */
    pRB getRBByName(const char *name) const {
        if (name == nullptr) return pRB();
        for (int i = (int)RBs.size() - 1; i >= 0; i--)
            if (strcmp(name, RBs[i]->name.c_str()) == 0) return RBs[i];
        return pRB();
    }

    /**
     * This method returns the shared pointer reference to the joint whose name
     * matches, or nullptr if it is not found
     */
    pRBJ getJointByName(const char *name) const {
        if (name == nullptr) return pRBJ();
        for (int i = (int)joints.size() - 1; i >= 0; i--)
            if (strcmp(name, joints[i]->name.c_str()) == 0) return joints[i];
        return pRBJ();
    }

    inline int getJointCount() { return (int)joints.size(); }

    inline int getRBCount() { return (int)RBs.size(); }

    int idx(const pRB& rb) const {
        if (rb.get() == nullptr)
            return -1;
        return rbIndices.at(rb);
    }

    int idx(const pRBJ& j) const {
        if (j.get() == nullptr)
            return -1;
        return jIndices.at(j);
    }

    //sets up data structures, indices and such, as needed...
    virtual void finalize() {
        rbIndices.clear();
        jIndices.clear();

        for (uint i = 0; i < RBs.size(); i++)
            rbIndices[RBs[i]] = i;

        for (uint i = 0; i < joints.size(); i++)
            jIndices[joints[i]] = i;

    }

    /**
     * returns a nullptr if no RBs are hit by the ray, otherwise returns the one that's
     * hit and is closest to the origin of the ray
     */
    pRB getFirstRBHitByRay(const Ray& ray){
        P3D p;
        return getFirstRBHitByRay(ray, p);
    }

    pRB getFirstRBHitByRay(const Ray &ray, P3D &intersectionPoint) {
        pRB selectedRB = pRB();
        double t = DBL_MAX;

        for (uint i = 0; i < RBs.size(); i++) {
            P3D p = RBs[i]->getRayHit(ray);
            if (p.isAtInfinity() == false) {
                double tTmp = ray.getRayParameterFor(p);
                if (tTmp < t) {
                    selectedRB = RBs[i];
                    t = tTmp;
                    intersectionPoint = p;
                }
            }
        }

        return selectedRB;
    }

    virtual void draw() {
        for (uint i = 0; i < RBs.size(); i++) {
            pRB rb = RBs[i];

            // draw collsion primitives
            if (df.showCollisionPrimitives)
                RBSRenderer::drawCollisionPrimitives(rb);

            // and end effectors
            if (df.showEndEffectors)
                RBSRenderer::drawEndEffectors(rb);

            RBSRenderer::draw3DModels(rb);
        }

        //draw moment of inertia last, as it is transparent
        if (df.showMOI)
            for (uint i = 0; i < RBs.size(); i++)
                RBSRenderer::drawMOI(RBs[i]);
    }

};
