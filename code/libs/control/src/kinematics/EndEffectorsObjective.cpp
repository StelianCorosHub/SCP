#include <control/kinematics/EndEffectorsObjective.h>
#include <optimization/utils.h>

EndEffectorsObjective::EndEffectorsObjective(KinematicMotionPlan *mp, double weight) {
    kmPlan = mp;
    this->w = weight;
}

double EndEffectorsObjective::evaluate() const {
    double retVal = 0.0;

    for (EndEffector &endEffector : kmPlan->endEffectors) {
        kmPlan->gcRR.setQ(kmPlan->robotPoses[endEffector.sIndex]);

        if (endEffector.constrainPosition()){
            P3D targetPos = endEffector.targetPosition;
            P3D currentPos = kmPlan->gcRR.p_of_q(endEffector.localCoordinates, endEffector.qIdxOfParent);
            V3D errPos(targetPos, currentPos);

            for (int k = 0; k < 3; k++)
                if (endEffector.positionMask[k] > 0)
                    retVal += 0.5 * errPos[k] * errPos[k] * endEffector.positionMask[k] * endEffector.posWeight * w;
        }

        if (endEffector.constrainOrientation()) {
            for (int i = 0; i < 3; i++) {
                if (endEffector.orientationMask[i] > 0) {
                    V3D targetVec = endEffector.targetFrame[i];
                    V3D currentVec = kmPlan->gcRR.v_of_q(endEffector.localFrame[i], endEffector.qIdxOfParent);
                    V3D errVec = targetVec - currentVec;
                    retVal += 0.5 * errVec.dot(errVec) * endEffector.orientationMask[i] * endEffector.rotWeight * w;
                }
            }
        }
    }

    return retVal;
}

void EndEffectorsObjective::addGradientTo(Vector& grad) const {
    Matrix J;

    // if root configuration is not part of the optimization, ignore the first 6 parameters
    int p_active_start_idx = (kmPlan->freezeRootDOFs == false) ? 0 : 6;

    for (EndEffector &endEffector : kmPlan->endEffectors) {
        kmPlan->gcRR.setQ(kmPlan->robotPoses[endEffector.sIndex]);

        if (endEffector.constrainPosition()) {
            P3D targetPos = endEffector.targetPosition;
            P3D currentPos = kmPlan->gcRR.p_of_q(endEffector.localCoordinates, endEffector.qIdxOfParent);
            V3D errPos(targetPos, currentPos);

            kmPlan->gcRR.compute_dpdq(endEffector.localCoordinates, endEffector.qIdxOfParent, J);

            //dO/dq = dO/dp * dp / dq, O(q) = 0.5 * v' * v, where v = p(q) - p_t
            for (int k = 0; k < 3; ++k) {
                if (endEffector.positionMask[k] > 0)
                    for (int l = p_active_start_idx; l < kmPlan->gcRR.getQSize(); l++) {
                        grad[endEffector.sIndex * kmPlan->getStateDim() + l - p_active_start_idx] +=
                                errPos[k] * J(k, l) * endEffector.positionMask[k] * endEffector.posWeight * w;
                }
            }
        }

        if (endEffector.constrainOrientation()) {
            for (int i = 0; i < 3; i++)
                if (endEffector.orientationMask[i] > 0) {
                    V3D targetVec = endEffector.targetFrame[i];
                    V3D currentVec = kmPlan->gcRR.v_of_q(endEffector.localFrame[i], endEffector.qIdxOfParent);
                    V3D errVec = targetVec - currentVec;
                    kmPlan->gcRR.compute_dvdq(endEffector.localFrame[i], endEffector.qIdxOfParent, J);
                    //dO/dq = dO/dp * dp / dq, O(q) = 0.5 * V' * V, where V = v(q) - v_t
                    for (int k = 0; k < 3; ++k) {
                        for (int l = p_active_start_idx; l < kmPlan->gcRR.getQSize(); l++) {
                            grad[endEffector.sIndex * kmPlan->getStateDim() + l - p_active_start_idx] +=
                                    -errVec[k] * J(k, l) * endEffector.orientationMask[i] * endEffector.rotWeight * w;
                        }
                    }
                }
        }
    }
}

void EndEffectorsObjective::addHessianEntriesTo(std::vector<SMTriplet> &hessianEntries) const {
    // if root configuration is not part of the optimization, ignore the first 6  parameters
    int p_active_start_idx = (kmPlan->freezeRootDOFs == false) ? 0 : 6;
    Matrix J;
    Matrix dJdql;

    for (EndEffector &endEffector : kmPlan->endEffectors) {
        kmPlan->gcRR.setQ(kmPlan->robotPoses[endEffector.sIndex]);

        if (endEffector.constrainPosition()) {
            P3D targetPos = endEffector.targetPosition;
            P3D currentPos = kmPlan->gcRR.p_of_q(endEffector.localCoordinates, endEffector.qIdxOfParent);
            V3D errPos(targetPos, currentPos);

            kmPlan->gcRR.compute_dpdq(endEffector.localCoordinates, endEffector.qIdxOfParent, J);

            for (int l = p_active_start_idx; l < kmPlan->gcRR.getQSize(); ++l) {
                bool hasNonZeros = false;
                if (kmPlan->use2ndOrderEETerms)
                    hasNonZeros = kmPlan->gcRR.compute_ddpdq_dqi(endEffector.localCoordinates,
                                                                 endEffector.qIdxOfParent,dJdql, l);
                for (int m = l; m < kmPlan->gcRR.getQSize(); m++) {
                    double val = 0.0;
                    for (int k = 0; k < 3; ++k) {
                        val += J(k, m) * J(k, l) * endEffector.positionMask[k] * endEffector.posWeight * w;

                        if (hasNonZeros && kmPlan->use2ndOrderEETerms)
                            val += errPos[k] * dJdql(k, m) *
                                   endEffector.positionMask[k]  * endEffector.posWeight * w;
                    }
                    ADD_HES_ELEMENT(hessianEntries,
                                    endEffector.sIndex * kmPlan->getStateDim() + l - p_active_start_idx,
                                    endEffector.sIndex * kmPlan->getStateDim() + m - p_active_start_idx,
                                    val);
                }
            }
        }

        if (endEffector.constrainOrientation()) {
            for (int i = 0; i < 3; ++i) {
                if (endEffector.orientationMask[i] > 0) {
                    // error
                    V3D targetVec = endEffector.targetFrame[i];
                    V3D currentVec = kmPlan->gcRR.v_of_q(endEffector.localFrame[i], endEffector.qIdxOfParent);
                    V3D errVec = targetVec - currentVec;
                    kmPlan->gcRR.compute_dvdq(endEffector.localFrame[i], endEffector.qIdxOfParent, J);

                    for (int l = p_active_start_idx; l < kmPlan->gcRR.getQSize(); ++l) {
                        bool hasNonZeros = false;
                        if (kmPlan->use2ndOrderEETerms)
                            hasNonZeros = kmPlan->gcRR.compute_ddvdq_dqi(endEffector.localFrame[i],
                                                                         endEffector.qIdxOfParent,dJdql, l);
                        for (int m = l; m < kmPlan->gcRR.getQSize(); ++m) {
                            double val = 0.0;
                            val += J.col(l).dot(J.col(m)) * endEffector.orientationMask[i] * endEffector.rotWeight * w;

                            if (hasNonZeros && kmPlan->use2ndOrderEETerms)
                                val += -dJdql.col(m).dot(errVec) * endEffector.orientationMask[i]  * endEffector.rotWeight * w;

                            ADD_HES_ELEMENT(hessianEntries,
                                            endEffector.sIndex * kmPlan->getStateDim() + l - p_active_start_idx,
                                            endEffector.sIndex * kmPlan->getStateDim() + m - p_active_start_idx,
                                            val);
                        }
                    }
                }
            }
        }
    }
}
