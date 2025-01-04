#include <RBSim/AFState.h>
#include <RBSim/ArticulatedFigure.h>
#include <RBSim/RBUtils.h>

AFState::AFState(const ArticulatedFigure& af){
    *this = af.getState();
}

AFState::AFState(ArticulatedFigure* af){
    *this = af->getState();
}

AFState::AFState(std::shared_ptr<ArticulatedFigure> pAF){
    *this = pAF->getState();
}

AFPose::AFPose(const ArticulatedFigure& af) {
    *this = af.getPose();
}

AFPose::AFPose(ArticulatedFigure* af) {
    *this = af->getPose();
}

AFPose::AFPose(std::shared_ptr<ArticulatedFigure> pAF) {
    *this = pAF->getPose();
}

void AFState::estimateQDotBasedOnDifferenceTo(const AFState& other, double dt){
    rootVel = V3D(rootPos, other.rootPos) / dt;
    rootAngVel = estimateAngularVelocity(rootQ, other.rootQ, dt);
    for (int i = 0; i < qDot.size(); i++)
        qDot[i] = (other.q[i] - q[i]) / dt;
}

