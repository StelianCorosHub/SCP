#include <control/kinematics/JointLimitsObjective.h>
#include <optimization/utils.h>

double JointLimitsObjective::evaluate() const {
    double retVal = 0;

    for (int i = 0; i < kmPlan->getNumberOfMotionSamples(); i++) {
        for (auto j : kmPlan->robot->joints) {
            if (j->hasJointLimits() == false)
                continue;
            int qIdx = kmPlan->robot->idx(j) + 6;

            double jointAngle = kmPlan->robotPoses[i][qIdx];

            double eps = fabs(j->maxAngle - j->minAngle) * 0.1;
            retVal += SoftBoundConstraint(j->minAngle - eps, j->maxAngle + eps, 1, 0.1)
                          .evaluate(jointAngle);
        }
    }

    return retVal * w;
}

void JointLimitsObjective::addGradientTo(Vector &grad) const {
    int p_active_start_idx = (kmPlan->freezeRootDOFs == false) ? 0 : 6;

    for (int i = 0; i < kmPlan->getNumberOfMotionSamples(); i++) {
        for (auto j : kmPlan->robot->joints) {
            if (j->hasJointLimits() == false)
                continue;
            int qIdx = kmPlan->robot->idx(j) + 6;

            double jointAngle = kmPlan->robotPoses[i][qIdx];

            double eps = fabs(j->maxAngle - j->minAngle) * 0.1;
            double dVal = SoftBoundConstraint(j->minAngle - eps, j->maxAngle + eps, 1, 0.1)
                              .computeDerivative(jointAngle);

            grad[kmPlan->getStateDim() * i + qIdx - p_active_start_idx] += dVal * w;
        }
    }
}

void JointLimitsObjective::addHessianEntriesTo(Array<SMTriplet> &hessianEntries) const {
    int p_active_start_idx = (kmPlan->freezeRootDOFs == false) ? 0 : 6;

    for (int i = 0; i < kmPlan->getNumberOfMotionSamples(); i++) {
        for (auto j : kmPlan->robot->joints) {
            if (j->hasJointLimits() == false)
                continue;
            int qIdx = kmPlan->robot->idx(j) + 6;

            double jointAngle = kmPlan->robotPoses[i][qIdx];

            double eps = fabs(j->maxAngle - j->minAngle) * 0.1;
            double ddVal = SoftBoundConstraint(j->minAngle - eps, j->maxAngle + eps, 1, 0.1)
                    .computeSecondDerivative(jointAngle);

            int hIdx = kmPlan->getStateDim() * i + qIdx - p_active_start_idx;

            ADD_HES_ELEMENT(hessianEntries, hIdx, hIdx, ddVal * w);
        }
    }
}
