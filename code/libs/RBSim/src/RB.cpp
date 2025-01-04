#include <RBSim/RB.h>
#include <utils/utils.h>

int tmpRB = 0;

#define INF_MASS 1e10

/**
 * Default constructor
 */
RB::RB(void) {
    Logger::print("creating an RB: %d\n", ++tmpRB);
}

/**
 * Default destructor
 */
RB::~RB(void) {
    Logger::print("deleting an RB: %d\n", --tmpRB);
}

/**
 *  returns the world coords moment of inertia of the rigid body.
 */
Matrix3x3 RB::MOI() {
    if (isStatic)
        return Matrix3x3::Identity() * INF_MASS;

    //R takes vectors from local coordinates to world coordinates
    Matrix3x3 R = orientation.toRotationMatrix();
    return R * MOI_local * R.transpose();
}

/**
 *  returns the inverse of the world coords moment of inertia of the rigid body.
 */
Matrix3x3 RB::invMOI(){
    if (isStatic)
        return Matrix3x3::Zero();
    //R takes vectors from local coordinates to world coordinates
    Matrix3x3 R = orientation.toRotationMatrix();
    return R * invMOI_local * R.transpose();
}

/**
 *  returns the local coords moment of inertia of the rigid body.
 */
Matrix3x3 RB::MOILocal(){
    return MOI_local;
}

/**
 *  returns the inverse of the local coords moment of inertia of the rigid body.
 */
Matrix3x3 RB::invMOILocal(){
    return invMOI_local;
}

double RB::invMass(){
    if (isStatic)
        return 0;
    else
        return 1.0 / mass;
}

/**
 * returns the 1 / mass of the rigid body
 **/
double RB::Mass(){
    if (isStatic)
        return INF_MASS;
    return mass;
}

P3D RB::getRayHit(const Ray &ray) {
    P3D intersectionPoint(DBL_MAX, DBL_MAX, DBL_MAX);
    double tMin = DBL_MAX;
    double t = tMin;

    RigidTransform RBTransform(orientation, position);

    // we will check all meshes
    for (uint i = 0; i < RBMeshes.size(); i++) {

        RigidTransform meshTransform = RBTransform * RBMeshes[i]->transform;

        RBMeshes[i]->model->position = meshTransform.T;
        RBMeshes[i]->model->orientation = meshTransform.R;

        if (RBMeshes[i]->model->hitByRay(ray.origin, ray.dir, t)) {
            if (t < tMin) {
                intersectionPoint = ray.origin + ray.dir * t;
                tMin = t;
            }
        }
    }

    return intersectionPoint;
}
