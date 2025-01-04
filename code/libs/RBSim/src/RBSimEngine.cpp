#include <RBSim/RBSimEngine.h>
#include <optimization/utils.h>

#include <RBSim/lcp/lcp.h>
#include <RBSim/lcp/memory.h>

void RBSimEngine::computeAndApplyConstraintForces(double h) {
    //gather up all constraints - from joints, motors, contacts/friction, etc.
    //then create a linear system A * lambda = b subject to inequality/lcp constraints,
    //and solve for lambda, the corresponding constraint forces

    constraintList.clear();

    for (auto j : joints)
        j->addUnboundedConstraintsToList(constraintList, h);

    //number of unbounded constraints
    int nu = constraintList.size();

    for (auto j : joints)
        j->addBoundedConstraintsToList(constraintList, h);

    //number of non-frictional contact constraints
    int nnfc = constraintList.size();
    int fnIdx = 0;

    //now we have the constraints arising from frictional contact
    collisionPoints = detectCollisions();

    for (auto& cp : collisionPoints)
        cp.addConstraintForNormalForce(constraintList, h);
    for (auto& cp : collisionPoints)
        cp.addConstraintsForFrictionForces(constraintList, h, nnfc + fnIdx++);

    //this is the number of constraints
    int n = (int)constraintList.size();
    if (n == 0) return;

    memoryManager.initialize(n);

    int nPad = memoryManager.nPad(n);
    double *A = memoryManager.A.get_ref();
    double *bx = memoryManager.bx.get_ref();
    double *lohi = memoryManager.lohi.get_ref();
    int *findex = memoryManager.findex.get_ref();

    //we will build the right hand side b of the system we need to solve to compute constraint forces:
    //b = -(CDot_t + h J*M^-1*F_ext + ERP / h * C)

    for (int i = 0; i < n; i++) {
        auto c = constraintList[i];
        //CDot_t = J * qDot
        bx[PBX__MAX * i + PBX_B]  = -(c.l_i.dot(c.rb_i->velocity) + c.a_i.dot(c.rb_i->angularVelocity)
                 + c.l_j.dot(c.rb_j->velocity) + c.a_j.dot(c.rb_j->angularVelocity));

        //h J*M^-1*F_ext
        bx[PBX__MAX * i + PBX_B] += -h * (c.l_i.dot(c.rb_i->invMass() * rbForces[c.rb_i]) + c.a_i.dot(c.rb_i->invMOI() * rbTorques[c.rb_i])
                      + c.l_j.dot(c.rb_j->invMass() * rbForces[c.rb_j]) + c.a_j.dot(c.rb_j->invMOI() * rbTorques[c.rb_j]));

        //ERP / h * c
        bx[PBX__MAX * i + PBX_B] += -c.erp / h * c.c;

        //feedforward term
        bx[PBX__MAX * i + PBX_B] += c.cDotFF;

        bx[PBX__MAX * i + PBX_X] = 0;
    }

    //this is the A matrix: A = h J * M^-1 * J' + diag_reg
    for (int i1 = 0; i1 < n; i1++) {
        auto c1 = constraintList[i1];
        //compute the J * M^-1 part first
        V3D l_i_m = V3D(c1.l_i * c1.rb_i->invMass());
        V3D a_i_I = V3D((c1.a_i.transpose() * c1.rb_i->invMOI()).transpose());
        V3D l_j_m = V3D(c1.l_j * c1.rb_j->invMass());
        V3D a_j_I = V3D((c1.a_j.transpose() * c1.rb_j->invMOI()).transpose());

        //and now multiply by everything by J' - boils down to a bunch of dot products
        for (int i2 = 0; i2 <= i1; i2++) {
            auto c2 = constraintList[i2];
            double val = 0;

            if (c1.rb_i.get() == c2.rb_i.get())
                val += l_i_m.dot(c2.l_i) + a_i_I.dot(c2.a_i);
            if (c1.rb_i.get() == c2.rb_j.get())
                val += l_i_m.dot(c2.l_j) + a_i_I.dot(c2.a_j);
            if (c1.rb_j.get() == c2.rb_i.get())
                val += l_j_m.dot(c2.l_i) + a_j_I.dot(c2.a_i);
            if (c1.rb_j.get() == c2.rb_j.get())
                val += l_j_m.dot(c2.l_j) + a_j_I.dot(c2.a_j);

            val *= h;

            //TODO: perhaps writing in the lower diagonal block only would suffice...
            A[i1 * nPad + i2] = val;
            if (i1 != i2)
                A[i2 * nPad + i1] = val;
        }

        A[i1 * nPad + i1] += c1.cfr;
    }

    //now take care of the constraints we have on the constraint forces...
    for (int i = 0; i < n; i++) {
        auto c = constraintList[i];
        lohi[PLH__MAX * i + PLH_LO] = c.cf_min;
        lohi[PLH__MAX * i + PLH_HI] = c.cf_max;
        findex[i] = c.pcindex;
    }

    //solve the lcp problem
    solveLCP(&memoryManager, n, nu);

    //finally, compute the world-coordinates constraint forces / torques: F_c = J' * lambda
    for (int i = 0; i < n; i++) {
        double lambda_i = bx[PBX__MAX * i + PBX_X];
        auto c = constraintList[i];
        if (c.cf_val != nullptr)
            *c.cf_val = lambda_i;
        rbForces[c.rb_i] += c.l_i * lambda_i;
        rbForces[c.rb_j] += c.l_j * lambda_i;
        rbTorques[c.rb_i] += c.a_i * lambda_i;
        rbTorques[c.rb_j] += c.a_j * lambda_i;
    }

    //and done :)
}

void RBSimEngine::step(double dt) {
    //add gravitational forces to all RBs
    for (auto rb : RBs) {
        if (rb->isStatic == false) {
            rbForces[rb] += rb->Mass() * UP * -9.8;
            // although an abuse of notation, add the w x Iw term to the net torque here so that we don't
            // have to keep carrying it with us in downstream computations
            rbTorques[rb] -= rb->angularVelocity.cross(V3D(rb->MOI() * rb->angularVelocity));
        }
    }

    //add feedforward motor torques, too...
    for (auto j : joints)
        j->applyFFMotorTorques(rbTorques);

    //compute constraint forces / torques
    computeAndApplyConstraintForces(dt);

    //and now integrate state forward in time using a symplectic euler time-stepping scheme
    for (auto rb : RBs) {
        V3D netForce = rbForces[rb];
        V3D netTorque = rbTorques[rb];

        rb->velocity = rb->velocity + dt * rb->invMass() * netForce;
        rb->position = rb->position + dt * rb->velocity;

        rb->angularVelocity = rb->angularVelocity + dt * rb->invMOI() * netTorque;
        rb->orientation = integrateOrientationForwardInTime(rb->orientation, rb->angularVelocity, dt);

/*      - below is the version that computes dq = 1/2 (0, w) * q_t and adds it to q_t
        Quaternion dq; dq.w() = 0; dq.vec() = 0.5 * rb->angularVelocity;
        dq *= rb->orientation;
        rb->orientation.w() += dq.w() * dt;
        rb->orientation.x() += dq.x() * dt;
        rb->orientation.y() += dq.y() * dt;
        rb->orientation.z() += dq.z() * dt;
        rb->orientation.normalize();
*/
        //finally, clear all forces and torques for this rb, since they've been used up
        rbForces[rb] = V3D();
        rbTorques[rb] = V3D();
    }
}

/**
 * We need an oracle that tells us if collision checking should be performed between rbA and rbB.
 * For now, we only check collisions if either rbA or rbB is static (but not both), and if there
 * is no joint connecting them.
 */

bool RBSimEngine::shouldCheckCollisionsBetween(pRB rbA, pRB rbB) {
    if (rbA == rbB)
        return false;

    if (checkCollisionsOnlyWithStaticObjects) {
        if (rbA->isStatic == false && rbB->isStatic == false)
            return false;

        if (rbA->isStatic && rbB->isStatic)
            return false;
    }

    for (auto j : joints) {
        if (j->parent == rbA && j->child == rbB)
            return false;
        if (j->parent == rbB && j->child == rbA)
            return false;
    }

    return true;
}


/**
 * perform collision detection, and return a list of all contact points between all pairs of rigid bodies
 */
CollisionPointList RBSimEngine::detectCollisions() {
    CollisionPointList res;

    //for now we check collisions only between static and non-static RBs
    for (uint i = 0; i < RBs.size(); i++) {
        for (uint j = i + 1; j < RBs.size(); j++) {
            if (shouldCheckCollisionsBetween(RBs[i], RBs[j]) == false)
                continue;
            for (auto cpA : RBs[i]->collisionShapes) {
                for (auto cpB: RBs[j]->collisionShapes) {
                    Array<ContactPoint> CPs = cpA->collideWith(cpB.get());
                    for (auto c: CPs)
                        res.push_back(RBCollisionObject(RBs[i], RBs[j], c));
                }
            }
        }
    }

    return res;
}


