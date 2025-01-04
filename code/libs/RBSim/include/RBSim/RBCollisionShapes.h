#pragma once

#include <RBSim/CollisionChecker.h>

class RB;

/**
 * Abstract class for collision primitives.
 */
class RBCollisionShape {
protected:
    RB* parentRB;
public:
    RBCollisionShape(RB* rb) {this->parentRB = rb;}
    virtual Array<ContactPoint> collideWith(RBCollisionShape* other) = 0;
    virtual void draw(const V3D& color, double alpha) = 0;
    virtual ~RBCollisionShape() {}
};

typedef std::shared_ptr<RBCollisionShape> pRBCollisionShape;

/**
 * Collision primitive shaped a sphere.
 */
class RBCollisionSphere : public RBCollisionShape {
public:
    //sphere in the local coordiante frame of the parent RB
    Sphere s;

    RBCollisionSphere(RB* rb, const P3D& pos, double radius) : RBCollisionShape(rb) {
        s.r = radius;
        s.p = pos;
    }

    ~RBCollisionSphere() {}

    Sphere getWorldCoordinatesPrimitive();
    virtual Array<ContactPoint> collideWith(RBCollisionShape* other);
    virtual void draw(const V3D& color, double alpha);

public:

};

/**
 * Collision primitive for an infinite plane
 */
class RBCollisionPlane : public RBCollisionShape {
public:
    //the plane, in the local coordinate frame of the parent RB
    Plane p;

    RBCollisionPlane(RB* rb, P3D p, V3D n) : RBCollisionShape(rb) {
        this->p.p = p;
        this->p.n = n;
    }

    ~RBCollisionPlane() {}

    Plane getWorldCoordinatesPrimitive();
    virtual Array<ContactPoint> collideWith(RBCollisionShape* other);
    virtual void draw(const V3D& color, double alpha);

public:

};

/**
 * Collision primitive for an infinite plane
 */
class RBCollisionBox : public RBCollisionShape {
public:
    //the box, expressed in the local coordinate frame of the parent RB
    Box3D b;

    RBCollisionBox(RB* rb, const P3D& p, const Quaternion& q, const V3D& dims) : RBCollisionShape(rb) {
        this->b.p = p;
        this->b.q = q;
        this->b.dims = dims;
    }

    ~RBCollisionBox() {}

    Box3D getWorldCoordinatesPrimitive();
    virtual Array<ContactPoint> collideWith(RBCollisionShape* other);
    virtual void draw(const V3D& color, double alpha);

public:

};
