#pragma once

#include <utils/mathUtils.h>

#include <utils/geoms.h>
#include <gui/model.h>
#include <RBSim/RBCollisionShapes.h>

/**
 * This class holds the 3d models we'll be using for visualization
 */
class RB3DModel {
public:

    RB3DModel(const std::string& fName){
        name = fName;
        if (!fileExists(name.c_str()))
            name = SCP_DATA_FOLDER "/" + fName;
        model = new Model(name);
    }

    ~RB3DModel() {
        delete model;
    }

    Model *model = nullptr;

    //transform used to position the 3d model in the local coordinate frame of the rigid body
    RigidTransform transform;
    std::string name;
    // color for the model, as well as its associated alpha value
    V3D color = V3D(0.5, 0.5, 0.5);
    double alpha = 1.0;
};

typedef std::shared_ptr<RB3DModel> pRB3DModel;

/**
 * End effector point. It's assumed to be a sphere with a certain radius.
 */
class RBEndEffector {
public:
    double radius = 0.01;
    // center of ee in local frame
    P3D pos = P3D(0, 0, 0);
    // with ground, obstacle etc.
    bool inContact = false;
};

/**
 * Rigid body class.
 */
class RB {
private:
    // mass
    double mass = 1.0;
    // moment of inertia of the rigid body, expressed in its local coordinate frame
    Matrix3x3 MOI_local = Matrix3x3::Identity() * 0.00166;
    //store the inverse moi in local coordinates for quick access...
    Matrix3x3 invMOI_local = Matrix3x3::Identity() * (1.0 / 0.00166);

public:
    std::string name;

/*
 * PHYSICAL PROPERTIES
 */

    // coefficient of restitution that determines how bouncy the contacts are
    double restitutionCoeff = 0;
    // friction coefficient
    double frictionCoeff = 0.8;

    // collision primitives that are relevant for this rigid body
    Array<pRBCollisionShape> collisionShapes;

    // end effector points
    Array<RBEndEffector> endEffectorPoints;

    // if this body is frozen in space, i.e. it is static, that's
    // equivalent to it having zero velocity and infinite mass / moment of inertia (no force will move it)
    bool isStatic = false;

/*
 * STATE VARIABLES
 */

    // the position of the center of mass of the rigid body, in world coords. In
    // local coordinates this corresponds to the point (0,0,0) aka origin.
    P3D position;
    // its orientation - rotates from the local coordinate frame of the rigid
    // body to world coordinates
    Quaternion orientation;
    // the linear velocity of the center of mass, in world coords
    V3D velocity = V3D(0, 0, 0);
    // the angular velocity in world coords
    V3D angularVelocity = V3D(0, 0, 0);

/*
 * VISUALIZATION PARAMETERS
 */

    // managed via mouse interactions in the GUI
    bool selected = false;

    // meshes used to visualize the rigid body, also used for picking
    Array<pRB3DModel> RBMeshes;

public:
    /**
     * Default constructor
     */
    RB(void);

    /**
     * Default destructor
     */
    virtual ~RB(void);

    /**
     * This method returns the coordinates of the point that is passed in as a
     * parameter(expressed in local coordinates), in world coordinates.
     */
    inline P3D getWorldCoordinates(const P3D &pLocal) const {
        // pWorld = pos + R * V3D(origin, pLocal)
        return position + orientation * V3D(P3D(), pLocal);
    }

    /**
     * This method returns the vector that is passed in as a parameter(expressed
     * in local coordinates), in world coordinates.
     */
    inline V3D getWorldCoordinates(const V3D &vLocal) const {
        // the rigid body's orientation is a unit quaternion. Using this, we can
        // obtain the global coordinates of a local vector
        return orientation * vLocal;
    }

    /**
     * This method is used to return the local coordinates of the point that is
     * passed in as a parameter (expressed in global coordinates)
     */
    inline P3D getLocalCoordinates(const P3D &pWorld) {
        return P3D() + orientation.inverse() * (V3D(position, pWorld));
    }

    /**
     * This method is used to return the local coordinates of the vector that is
     * passed in as a parameter (expressed in global coordinates)
     */
    inline V3D getLocalCoordinates(const V3D &vWorld) {
        // the rigid body's orientation is a unit quaternion. Using this, we can
        // obtain the global coordinates of a local vector
        return orientation.inverse() * vWorld;
    }

    /**
     * This method returns the world coordinates velocity for a p that is
     * expressed in the local coordinate frame of the rigid body. Evidently,
     * this point has 0 velocity in the local coordinate frame, as the rigid
     * body does not change its shape over time.
     */
    inline V3D getVelocityForLocalCoordsPoint(const P3D &pLocal) {
        // we need to compute the vector r, from the origin of the body to the
        // point of interest
        V3D r = V3D(pLocal);
        // the velocity is given by w x r + v. w and v are already
        // expressed in world coordinates, so we need to express r in world
        // coordinates too.
        return angularVelocity.cross(getWorldCoordinates(r)) + velocity;
    }

    /**
     *  returns the world coords moment of inertia of the rigid body.
     */
    Matrix3x3 MOI();

    /**
     *  returns the inverse of the world coords moment of inertia of the rigid body.
     */
    Matrix3x3 invMOI();

    /**
     *  returns the local coords moment of inertia of the rigid body.
     */
    Matrix3x3 MOILocal();

    /**
     *  returns the inverse of the local coords moment of inertia of the rigid body.
     */
    Matrix3x3 invMOILocal();

    /**
     * returns the 1 / mass of the rigid body
     * */
    double invMass();

    /**
     * returns the 1 / mass of the rigid body
     **/
    double Mass();

    void setMass(double m){
        this->mass = m;
    }

    /**
     * set the moment of inertia of the rigid body - symmetric 3x3 matrix, so
     * we need the six values for it.
     */
    inline void setMOI(double moi00, double moi11, double moi22, double moi01,
                       double moi02, double moi12) {
        MOI_local << moi00, moi01, moi02, moi01, moi11, moi12, moi02, moi12,
                moi22;
        invMOI_local = MOI_local.inverse();
    }
/*
    Matrix3x3 getMOI(const Matrix3x3& RToWorld) {
        return RToWorld * MOI_local * RToWorld.transpose();
    }

    Matrix3x3 getMOI(const Quaternion& q) {
        return getMOI(q.toRotationMatrix());
    }
*/
    /**
     *  returns the first point on the rigid body that is hit by the ray,
     *  or a point of infinity if the ray does not hit the rigid body.
     */
    P3D getRayHit(const Ray &ray);
};

typedef std::shared_ptr<RB> pRB;

inline pRB createBoxRB(const std::string name = "box", double lx = 0.1, double ly = 0.1, double lz = 0.1, double density = 1000) {
    pRB rb = pRB(new RB());
    rb->name = name;
    rb->setMass(lx * ly * lz * density);
    rb->setMOI(1.0 / 12 * rb->Mass() * (ly * ly + lz * lz), 1.0 / 12 * rb->Mass() * (lx * lx + lz * lz), 1.0 / 12 * rb->Mass() * (lx * lx + ly * ly), 0, 0, 0);

    rb->RBMeshes.push_back(std::make_shared<RB3DModel>("/meshes/cube.obj"));
    rb->RBMeshes.back()->model->scale = V3D(lx, ly, lz);

/*
    rb->collisionShapes.push_back(std::make_shared<RBCollisionSphere>(rb.get(), P3D(lx / 2.0, ly / 2.0, lz / 2.0), 0.01));
    rb->collisionShapes.push_back(std::make_shared<RBCollisionSphere>(rb.get(), P3D(lx / 2.0, -ly / 2.0, lz / 2.0), 0.01));
    rb->collisionShapes.push_back(std::make_shared<RBCollisionSphere>(rb.get(), P3D(-lx / 2.0, -ly / 2.0, lz / 2.0), 0.01));
    rb->collisionShapes.push_back(std::make_shared<RBCollisionSphere>(rb.get(), P3D(-lx / 2.0, ly / 2.0, lz / 2.0), 0.01));

    rb->collisionShapes.push_back(std::make_shared<RBCollisionSphere>(rb.get(), P3D(lx / 2.0, ly / 2.0, -lz / 2.0), 0.01));
    rb->collisionShapes.push_back(std::make_shared<RBCollisionSphere>(rb.get(), P3D(lx / 2.0, -ly / 2.0, -lz / 2.0), 0.01));
    rb->collisionShapes.push_back(std::make_shared<RBCollisionSphere>(rb.get(), P3D(-lx / 2.0, -ly / 2.0, -lz / 2.0), 0.01));
    rb->collisionShapes.push_back(std::make_shared<RBCollisionSphere>(rb.get(), P3D(-lx / 2.0, ly / 2.0, -lz / 2.0), 0.01));
*/

    rb->collisionShapes.push_back(std::make_shared<RBCollisionBox>(rb.get(), P3D(), Quaternion(), V3D(lx, ly, lz)));

    return rb;
}

inline pRB createStaticPlaneRB(const std::string name = "groundPlane", const P3D& p = P3D(0,0,0), const V3D& n = V3D(0,1,0)) {
    pRB rb = pRB(new RB());
    rb->name = name;
    rb->isStatic = true;

    rb->collisionShapes.push_back(std::make_shared<RBCollisionPlane>(rb.get(), p, n));

    return rb;
}

inline pRB createSphereRB(const std::string name = "sphere", double r = 0.1, double density = 1000) {
    pRB rb = pRB(new RB());
    rb->name = name;
    double m = 4.0 / 3 * r * r * r * M_PI * density;
    double I = 2.0 / 5 * m * r * r;
    rb->setMass(m);
    rb->setMOI(I, I, I, 0, 0, 0);

    rb->RBMeshes.push_back(std::make_shared<RB3DModel>("/meshes/sphere.obj"));
    rb->RBMeshes.back()->model->scale = V3D(2 * r, 2 * r, 2 * r);

    rb->collisionShapes.push_back(std::make_shared<RBCollisionSphere>(rb.get(), P3D(), r));

    return rb;
}

