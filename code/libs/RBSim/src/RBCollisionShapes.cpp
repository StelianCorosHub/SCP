#include <RBSim/RBCollisionShapes.h>
#include <RBSim/RB.h>
#include <gui/renderer.h>

Sphere RBCollisionSphere::getWorldCoordinatesPrimitive() {
    return Sphere(parentRB->getWorldCoordinates(s.p), s.r);
}

void RBCollisionSphere::draw(const V3D& color, double alpha){
    Sphere sW = getWorldCoordinatesPrimitive();
    drawSphere(sW.p, sW.r, color, alpha);
}

Array<ContactPoint> RBCollisionSphere::collideWith(RBCollisionShape* other){
    Sphere s = getWorldCoordinatesPrimitive();
    if (auto sOther = dynamic_cast<RBCollisionSphere*>(other))
        return collideSphereWithSphere(s, sOther->getWorldCoordinatesPrimitive());
    if (auto bOther = dynamic_cast<RBCollisionBox*>(other))
        return collideSphereWithBox(s, bOther->getWorldCoordinatesPrimitive());
    if (auto pOther = dynamic_cast<RBCollisionPlane*>(other))
        return collideSphereWithPlane(s, pOther->getWorldCoordinatesPrimitive());
    return Array<ContactPoint>();
}

Plane RBCollisionPlane::getWorldCoordinatesPrimitive() {
    return Plane(parentRB->getWorldCoordinates(p.p), parentRB->getWorldCoordinates(p.n));
}

void RBCollisionPlane::draw(const V3D& color, double alpha){
    Plane pW = getWorldCoordinatesPrimitive();
    drawRectangle(pW.p, pW.n, 0, Vector2d(100, 100), color, alpha);
}

Array<ContactPoint> RBCollisionPlane::collideWith(RBCollisionShape* other){
    Plane p = getWorldCoordinatesPrimitive();
    if (auto sOther = dynamic_cast<RBCollisionSphere*>(other))
        return collidePlaneWithSphere(p, sOther->getWorldCoordinatesPrimitive());
    if (auto bOther = dynamic_cast<RBCollisionBox*>(other))
        return collidePlaneWithBox(p, bOther->getWorldCoordinatesPrimitive());
    return Array<ContactPoint>();
}

Box3D RBCollisionBox::getWorldCoordinatesPrimitive() {
    return Box3D(parentRB->getWorldCoordinates(b.p), parentRB->orientation * this->b.q, b.dims);
}

void RBCollisionBox::draw(const V3D& color, double alpha){
    Box3D bW = getWorldCoordinatesPrimitive();
    drawCuboid(bW.p, bW.q, bW.dims, color, alpha);
}

Array<ContactPoint> RBCollisionBox::collideWith(RBCollisionShape* other){
    Box3D b = getWorldCoordinatesPrimitive();
    if (auto sOther = dynamic_cast<RBCollisionSphere*>(other))
        return collideBoxWithSphere(b, sOther->getWorldCoordinatesPrimitive());
    if (auto bOther = dynamic_cast<RBCollisionBox*>(other))
        return collideBoxWithBox(b, bOther->getWorldCoordinatesPrimitive());
    if (auto pOther = dynamic_cast<RBCollisionPlane*>(other))
        return collideBoxWithPlane(b, pOther->getWorldCoordinatesPrimitive());
    return Array<ContactPoint>();
}


