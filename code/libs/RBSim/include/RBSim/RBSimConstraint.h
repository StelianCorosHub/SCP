#pragma once

#include <RBSim/RB.h>

#include <utils/mathUtils.h>
#include <utils/utils.h>

//class used to model a single scalar constraint in RBSimEngine
class RBSimConstraint {
public:
    //constraints are used to couple pairs of rigid bodies to each other. This constraint couples rb_i and rb_j to each other.
    pRB rb_i = pRB(), rb_j = pRB();

    //these are the components of the jacobian of the constraint wrt to the linear and angular DOFs of rigid bodies i and j
    Vector3d l_i, a_i, l_j, a_j;

    //the current value of the constraint, as well as the current time derivative of the constraint's value
    double c = 0;

    //this is a feedforward term used to set non-zero (feedback terms notwhistanding) targets for CDot after the timestep
    double cDotFF = 0;

    //this is the error reduction parameter, a scalar between 0 and 1, which indicates how much of the value of the constraint
    //we want to correct with each time step
    double erp = 0;

    //this is the constraint force regularizer parameter
    double cfr = 0;

    //these are min and max values for the constraint force.
    double cf_min = -INFINITY;
    double cf_max = INFINITY;

    //once a constraint force is computed, its value will be written into here (unless this pointer is null).
    double* cf_val = nullptr;

    //sometimes, the bounds for one constraint depend on the values of the force for another constraint (e.g. friction
    //forces are bounded by the value of the force in the normal direction). To deal with such cases, we will
    //store the index of the parent constraint in here.
    int pcindex = -1;

public:
    RBSimConstraint(pRB rb_i, pRB rb_j){
        this->rb_i = rb_i;
        this->rb_j = rb_j;
    }

};

