#pragma once

#include <RBSim/RB.h>
#include <RBSim/CollisionChecker.h>

#include <utils/mathUtils.h>
#include <utils/utils.h>

class RBCollisionObject {
public:
    pRB rbA, rbB;
    //the contact point, in world coordinates
    P3D p;
    //normal, also in world coordinates
    V3D n;
    //and the penetration depth
    double d;

    double kp = 2.4 * 1e6;
    double kd = 0.8 * 1e5;

    double fn = 0;
    double ft[2] = {0, 0};

    RBCollisionObject(pRB rbA, pRB rbB, const ContactPoint& cp) {
        this->rbA = rbA;
        this->rbB = rbB;
        this->p = cp.p;
        this->n = cp.n;
        this->d = cp.d;
    }

    double getRestitutionCoefficient() {
        return rbA->restitutionCoeff * 0.5 + rbB->restitutionCoeff * 0.5;
    }

    double getFrictionCoefficient() {
        return rbA->frictionCoeff * 0.5 + rbB->frictionCoeff * 0.5;
    }

    void addConstraintForNormalForce(Array<RBSimConstraint>& constraints, double sim_dt) {
        //we want to have zero relative velocity at the contact point, but in the normal direction the constraint
        //force needs to be positive (no pulling)

        V3D va = rbA->getVelocityForLocalCoordsPoint(rbA->getLocalCoordinates(p));
        V3D vb = rbB->getVelocityForLocalCoordsPoint(rbB->getLocalCoordinates(p));
        V3D vRel = va - vb;

        V3D ra(rbA->position, p);
        V3D rb(rbB->position, p);

        constraints.push_back(RBSimConstraint(rbA, rbB));
        constraints.back().l_i = -n; constraints.back().a_i = n.cross(ra);
        constraints.back().l_j = n; constraints.back().a_j = -n.cross(rb);

        constraints.back().c = d;
        constraints.back().cf_val = &fn;

        // ERP = h*kp / (h*kp + kd)
        constraints.back().erp = sim_dt * kp / (sim_dt * kp + kd);
        // CFM = 1 / (h*kp + kd)
        constraints.back().cfr = 1.0 / (sim_dt * kp + kd);

        //only allow pushing forces
        constraints.back().cf_min = 0;

        //add restitution...
        double cDot = vRel.dot(n);
        if (cDot > 0)
            constraints.back().cDotFF = -cDot * getRestitutionCoefficient();

    }

    void addConstraintsForFrictionForces(Array<RBSimConstraint>& constraints, double sim_dt, int ncIdx) {
        double mu = getFrictionCoefficient();

        if (mu <= 0)
            return;

        V3D ra(rbA->position, p);
        V3D rb(rbB->position, p);

        V3D tDir[2];
        getOrthonormalVectors(n, tDir[0], tDir[1]);

        for (int i = 0; i < 2; i++) {
            constraints.push_back(RBSimConstraint(rbA, rbB));
            constraints.back().l_i = tDir[i]; constraints.back().a_i = ra.cross(tDir[i]);
            constraints.back().l_j = -tDir[i]; constraints.back().a_j = -rb.cross(tDir[i]);

            constraints.back().cf_val = &ft[i];
            constraints.back().cfr = 1e-6;
            constraints.back().cf_min = -mu;
            constraints.back().cf_max = mu;
            constraints.back().pcindex = ncIdx;
        }
    }

};

typedef Array<RBCollisionObject> CollisionPointList;

