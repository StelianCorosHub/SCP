#include <control/kinematics/RegularizerObjective.h>
#include <optimization/utils.h>

PoseRegularizer::PoseRegularizer(KinematicMotionPlan *mp, double weight) {
    kmPlan = mp;
    this->w = weight;
}

double PoseRegularizer::evaluate() const {
    double retVal = 0;

    int p_active_start_idx = (kmPlan->freezeRootDOFs == false) ? 0 : 6;

    for (auto t : kmPlan->targetPoses)
        for (int i = p_active_start_idx; i < kmPlan->gcRR.getQSize(); i++) {
            double tmpV = (t.target[i] - kmPlan->robotPoses[t.sIndex][i]);
            retVal += 0.5 * tmpV * tmpV * t.jw[i] * w;
        }

    return retVal;
}

void PoseRegularizer::addGradientTo(Vector &grad) const {
    int p_active_start_idx = (kmPlan->freezeRootDOFs == false) ? 0 : 6;

    for (auto t : kmPlan->targetPoses)
        for (int i = p_active_start_idx; i < kmPlan->gcRR.getQSize(); i++) {
            double tmpV = (t.target[i] - kmPlan->robotPoses[t.sIndex][i]);
            grad[kmPlan->getStateDim() * t.sIndex + i - p_active_start_idx] += -tmpV * t.jw[i] * w;
        }
}

void PoseRegularizer::addHessianEntriesTo(Array<SMTriplet> &hessianEntries) const {
    int p_active_start_idx = (kmPlan->freezeRootDOFs == false) ? 0 : 6;

    for (auto t : kmPlan->targetPoses)
        for (int i = p_active_start_idx; i < kmPlan->gcRR.getQSize(); i++) {
            int idx = kmPlan->getStateDim() * t.sIndex + i - p_active_start_idx;
            ADD_HES_ELEMENT(hessianEntries, idx, idx, t.jw[i] * w);
        }
}

/*
SmoothMotionRegularizer::SmoothMotionRegularizer(KinematicMotionPlan *mp,
                                                 int startQIndex, int endQIndex,
                                                 std::string description,
                                                 double weight)
    : SubObjective(description, weight) {
    kmPlan = mp;
    this->startQIndex = startQIndex;
    this->endQIndex = endQIndex;
}

//////////////////////////////////////////////////////////////////////////////

double SmoothMotionRegularizer::evaluate(const dVector &x) const {
    double retVal = 0;

    for (int j = 1; j < kmPlan->nStateSamples - 1; j++) {
        for (int i = startQIndex; i <= endQIndex; i++) {
            // j means samples in time, i is the joint index...
            double a = kmPlan->robotStates[j - 1][i] -
                       2 * kmPlan->robotStates[j][i] +
                       kmPlan->robotStates[j + 1][i];
            retVal += 0.5 * a * a * weight;
        }
    }

    return retVal;
}

void SmoothMotionRegularizer::addGradientTo(const dVector &x,
                                            dVector &grad) const {
    int stateDim = kmPlan->getStateDim();
    int offset = (kmPlan->optimizeRootConfiguration) ? 0 : 6;
    for (int j = 1; j < kmPlan->nStateSamples - 1; j++) {
        for (int i = startQIndex; i <= endQIndex; i++) {
            double a = kmPlan->robotStates[j - 1][i] -
                       2 * kmPlan->robotStates[j][i] +
                       kmPlan->robotStates[j + 1][i];
            if (i >= offset) {
                grad[stateDim * (j - 1) + i - offset] += a * weight;
                grad[stateDim * (j + 0) + i - offset] += -2 * a * weight;
                grad[stateDim * (j + 1) + i - offset] += a * weight;
            }
        }
    }
}

void SmoothMotionRegularizer::addHessianEntriesTo(
    const dVector &x, std::vector<MTriplet> &hessianEntries) const {
    int stateDim = kmPlan->getStateDim();
    int offset = (kmPlan->optimizeRootConfiguration) ? 0 : 6;

    for (int j = 1; j < kmPlan->nStateSamples - 1; j++) {
        for (int i = startQIndex; i <= endQIndex; i++)
            if (i >= offset)
                for (int k1 = -1; k1 <= 1; k1++)
                    for (int k2 = k1; k2 <= 1; k2++) {
                        int val = 1;
                        if (k2 == 0 && k1 == 0)
                            val = 4;
                        else if (k1 == 0 || k2 == 0)
                            val = -2;
                        ADD_HES_ELEMENT(
                            hessianEntries, stateDim * (j + k1) + i - offset,
                            stateDim * (j + k2) + i - offset, val * weight);
                    }
    }
}

/////////////////////////////////////////////////////////////////////////

TargetStatesObjective::TargetStatesObjective(KinematicMotionPlan *mp,
                                             std::string description,
                                             double weight)
    : SubObjective(description, weight) {
    kmPlan = mp;
}

double TargetStatesObjective::evaluate(const dVector &x) const {
    int stateDim = kmPlan->getStateDim();
    double retVal = 0;

    for (uint j = 0; j < kmPlan->sTargets.size(); j++) {
        for (int i = 0; i < kmPlan->defaultPose.size(); i++) {
            double tmpV = (kmPlan->sTargets[j].target[i] -
                           kmPlan->robotStates[kmPlan->sTargets[j].sIndex][i]);
            retVal += 0.5 * tmpV * tmpV * weight * kmPlan->sTargets[j].weight;
        }
    }

    return retVal;
}

void TargetStatesObjective::addGradientTo(const dVector &x,
                                          dVector &grad) const {
    int stateDim = kmPlan->getStateDim();

    for (uint j = 0; j < kmPlan->sTargets.size(); j++) {
        for (int i = 0; i < kmPlan->defaultPose.size(); i++) {
            double tmpV = (kmPlan->sTargets[j].target[i] -
                           kmPlan->robotStates[kmPlan->sTargets[j].sIndex][i]);
            if (kmPlan->optimizeRootConfiguration)
                grad[stateDim * kmPlan->sTargets[j].sIndex + i] +=
                    -tmpV * weight * kmPlan->sTargets[j].weight;
            else if (i >= 6)
                grad[stateDim * kmPlan->sTargets[j].sIndex + i - 6] +=
                    -tmpV * weight * kmPlan->sTargets[j].weight;
        }
    }
}

void TargetStatesObjective::addHessianEntriesTo(
    const dVector &x, std::vector<MTriplet> &hessianEntries) const {
    int stateDim = kmPlan->getStateDim();

    for (uint j = 0; j < kmPlan->sTargets.size(); j++) {
        for (int i = 0; i < kmPlan->defaultPose.size(); i++) {
            if (kmPlan->optimizeRootConfiguration)
                ADD_HES_ELEMENT(hessianEntries,
                                stateDim * kmPlan->sTargets[j].sIndex + i,
                                stateDim * kmPlan->sTargets[j].sIndex + i,
                                weight * kmPlan->sTargets[j].weight);
            else if (i >= 6)
                ADD_HES_ELEMENT(hessianEntries,
                                stateDim * kmPlan->sTargets[j].sIndex + i - 6,
                                stateDim * kmPlan->sTargets[j].sIndex + i - 6,
                                weight * kmPlan->sTargets[j].weight);
        }
    }
}
*/

