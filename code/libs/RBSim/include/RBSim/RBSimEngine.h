#pragma once

#include <RBSim/RBSRepo.h>
#include <RBSim/RBJoint.h>
#include <RBSim/RB.h>
#include <RBSim/ArticulatedFigure.h>
#include <RBSim/RBSLoader.h>
#include <RBSim/RBSimConstraint.h>

#include <RBSim/CollisionChecker.h>
#include <RBSim/RBCollisionObject.h>

#include <utils/mathUtils.h>
#include <utils/utils.h>
#include <RBSim/lcp/lcp.h>

#include <gui/renderer.h>

/**
 * This class, which is a container for robots, rigid bodies and joints, is
 * a simulation world governed by Newton's laws of motion.
 */
class RBSimEngine : public RBSRepo {
private:
    std::unordered_map<pRB, V3D> rbForces;
    std::unordered_map<pRB, V3D> rbTorques;

    Array<V3D> constraintForces;
    Array<V3D> constraintTorques;

    Array<RBSimConstraint> constraintList;

    LCPMemoryManager memoryManager;

public:
    bool checkCollisionsOnlyWithStaticObjects = true;
    CollisionPointList collisionPoints;

public:
    // the constructor
    RBSimEngine() {}

    // the destructor
    virtual ~RBSimEngine() {

    }

    void addRB(pRB rb) {
        RBs.push_back(rb);
        rbIndices[RBs.back()] = RBs.size() - 1;
    }

    void addJoint(pRBJ j) {
        joints.push_back(j);
        jIndices[joints.back()] = joints.size() - 1;
    }

    void addRobot(pRobot robot) {
        for (auto rb : robot->RBs)
            addRB(rb);
        for (auto j : robot->joints)
            addJoint(j);
    }

    //this method adds a force f to the rigid body rb. The force is applied at pLocal,
    //where, as the name suggests, pLocal is expressed in the local coordinate frame
    //of rb.
    void addForceToRB(pRB rb, const V3D& f, const P3D& pLocal = P3D(0,0,0)) {
        rbForces[rb] += f;
        rbTorques[rb] += rb->getWorldCoordinates(V3D(P3D(), pLocal)).cross(f);
    }

    virtual void loadFromFile(const char *fName, bool loadVisuals = true) {
        // load robot from rbLoader
        RBSLoader rbsLoader(fName, loadVisuals);
        rbsLoader.populate(this);
    }

    /**
     * Advance the state of the simulation forward in time.
     */
    virtual void step(double dt);

    /**
     * perform collision detection, and return a list of all contact points between all pairs of rigid bodies
     */
    CollisionPointList detectCollisions();

    /**
     * We need an oracle that tells us if collision checking should be performed between rbA and rbB.
     * For now, we only check collisions if either rbA or rbB is static (but not both), and if there
     * is no joint connecting them.
     */
    bool shouldCheckCollisionsBetween(pRB rbA, pRB rbB);

    /**
     * Computes all constraint forces and applies them to the relevant RBs
     */
    void computeAndApplyConstraintForces(double dt);

    Array<double> getFullWorldState() {
        Array<double> s;
        for (auto rb : RBs){
            s.push_back(rb->position.x); s.push_back(rb->position.y); s.push_back(rb->position.z);
            s.push_back(rb->velocity.x()); s.push_back(rb->velocity.y()); s.push_back(rb->velocity.z());
            s.push_back(rb->orientation.w()); s.push_back(rb->orientation.x()); s.push_back(rb->orientation.y()); s.push_back(rb->orientation.z());
            s.push_back(rb->angularVelocity.x()); s.push_back(rb->angularVelocity.y()); s.push_back(rb->angularVelocity.z());
        }
        return s;
    }

    void setFullWorldState(const Array<double>& ws) {
        int idx = 0;
        for (auto rb : RBs) {
            rb->position.x = ws[idx++]; rb->position.y = ws[idx++]; rb->position.z = ws[idx++];
            rb->velocity.x() = ws[idx++]; rb->velocity.y() = ws[idx++]; rb->velocity.z() = ws[idx++];
            rb->orientation.w() = ws[idx++]; rb->orientation.x() = ws[idx++]; rb->orientation.y() = ws[idx++]; rb->orientation.z() = ws[idx++];
            rb->angularVelocity.x() = ws[idx++]; rb->angularVelocity.y() = ws[idx++]; rb->angularVelocity.z() = ws[idx++];
        }
        collisionPoints.clear();
    }

    /**
     * Draw rigid bodies and robots belong int the world.
     */
    inline void drawRBs() {
        //now we draw the individual RBs, in two passes, to account a little bit for
        for (uint i = 0; i < RBs.size(); i++){
            pRB rb = RBs[i];

            // Then draw collsion spheres
            if (df.showCollisionPrimitives)
                RBSRenderer::drawCollisionPrimitives(rb);

            // and end effectors
            if (df.showEndEffectors)
                RBSRenderer::drawEndEffectors(rb);
        }

        for (uint i = 0; i < RBs.size(); i++){
            pRB rb = RBs[i];

            // Then draw meshes (because of blending)
            if (df.showMeshes && df.showMOI == false) {
                if (df.showCollisionPrimitives || df.showEndEffectors)
                    RBSRenderer::draw3DModels(rb, 0.2);
                else
                    RBSRenderer::draw3DModels(rb);
            }

            // and now MOIs
            if (df.showMOI)
                RBSRenderer::drawMOI(rb);
        }
    }

    inline void drawContactForces(double scaleFactor = 0.01) {
        for (auto cp : collisionPoints) {
            V3D t1, t2;
            getOrthonormalVectors(cp.n, t1, t2);

            V3D f = cp.n * cp.fn + t1 * cp.ft[0] + t2 * cp.ft[1];
            drawArrow3d(cp.p, f * scaleFactor, 0.01);
        }
    }

};

