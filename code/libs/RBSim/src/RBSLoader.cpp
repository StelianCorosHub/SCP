#include "RBSim/RBSLoader.h"

RBSLoader::RBSLoader(const char *fName, bool loadGUIElements)
    : loadGUIElements(loadGUIElements) {

    this->verbose = (loadGUIElements == true);
    loadBufferFromFile(fName);
}

RBSLoader::~RBSLoader() {
    freeBuffer();
}

int RBSLoader::getRBSLineType() {
    int idx = 0;

    goToNextNonWhiteSpaceCharacter();

    if (endOfTextStream())
        return RBS_END_OF_FILE;

    if (peek() == '#')
        return RBS_COMMENT;

    while (idx < RBSKeywords.size()){
        if (parseIfKeyword(RBSKeywords[idx].keywordText.c_str()))
            return RBSKeywords[idx].keywordID;
        idx++;
    }
    return RBS_UNKNOWN;
}

void RBSLoader::populate(RBSRepo* rbsRepo){
    RB* rb = nullptr;
    RBJoint* j = nullptr;

    // this is where it happens.
    while (!endOfTextStream()) {

        int tokenID = getRBSLineType();

        switch (tokenID) {
            case RBS_RB:
                // create a new rigid body and load it up
                rbsRepo->RBs.push_back(pRB(createAndPopulateRB()));
                break;
            case RBS_JOINT:
                // create a new joint and load it up
                rbsRepo->joints.push_back(pRBJ(createAndPopulateRBJoint(rbsRepo)));
                break;
            case RBS_COMMENT:
                parseEntireLine();
                break;
            case RBS_UNKNOWN:
                parserError("Unknown token");
                parseEntireLine();
                break;
            case RBS_END_OF_FILE:
                break;
            default:
                parserError("Unexpected token %s", RBSKeywords[tokenID].keywordText.c_str());
                parseEntireLine();

        }
    }

    rbsRepo->finalize();
}

RB* RBSLoader::createAndPopulateRB() {
    RB* rb = new RB();

    while (!endOfTextStream()) {
        int tokenID = getRBSLineType();

        switch (tokenID) {
            case RBS_NAME:
                parseName(rb->name);
            break;
            case RBS_MESH_NAME:
                if (loadGUIElements) {
                    std::string meshName;
                    parseName(meshName);
                    rb->RBMeshes.push_back(std::make_shared<RB3DModel>(meshName));
                }
                break;
            case RBS_MESH_COLOR:
                if (loadGUIElements) {
                    double r = 0, g = 0, b = 0;
                    parseDouble(r);
                    parseDouble(g);
                    parseDouble(b);
                    rb->RBMeshes.back()->color = V3D(r, g, b);
                }
                break;
            case RBS_MESH_TRANSFORMATION:
                if (loadGUIElements) {
                    Quaternion q;
                    parseDouble(q.w()); parseDouble(q.x()); parseDouble(q.y()); parseDouble(q.z());
                    P3D p;
                    parseDouble(p.x); parseDouble(p.y); parseDouble(p.z);
                    if (rb->RBMeshes.size() > 0) {
                        rb->RBMeshes.back()->transform.R = q;
                        rb->RBMeshes.back()->transform.T = p;
                    }
                }
                break;
            case RBS_MASS: {
                    double m = 0;
                    parseDouble(m);
                    rb->setMass(m);
                }
                break;
            case RBS_MOI: {
                    double t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0, t6 = 0;
                    parseDouble(t1); parseDouble(t2); parseDouble(t3);
                    parseDouble(t4); parseDouble(t5); parseDouble(t6);
                    rb->setMOI(t1, t2, t3, t4, t5, t6);
                }
                break;
            case RBS_FRICTION_COEFF:
                parseDouble(rb->frictionCoeff);
                break;
            case RBS_REST_COEFF:
                parseDouble(rb->restitutionCoeff);
                break;
            case RBS_IS_STATIC:
                rb->isStatic = true;
                break;

            case RBS_POSITION:
                parseDouble(rb->position.x);
                parseDouble(rb->position.y);
                parseDouble(rb->position.z);
                break;
            case RBS_VELOCITY:
                parseDouble(rb->velocity.x());
                parseDouble(rb->velocity.y());
                parseDouble(rb->velocity.z());
                break;
            case RBS_ORIENTATION:
                parseDouble(rb->orientation.w());
                parseDouble(rb->orientation.x());
                parseDouble(rb->orientation.y());
                parseDouble(rb->orientation.z());
                break;
            case RBS_ANGULAR_VELOCITY:
                parseDouble(rb->angularVelocity.x());
                parseDouble(rb->angularVelocity.y());
                parseDouble(rb->angularVelocity.z());
                break;
            case RBS_COLLISION_BOX: {
                    P3D p;
                    Quaternion q;
                    V3D dims;
                    parseDouble(p.x); parseDouble(p.y); parseDouble(p.z);
                    parseDouble(q.w()); parseDouble(q.x()); parseDouble(q.y()); parseDouble(q.z());
                    parseDouble(dims.x()); parseDouble(dims.y()); parseDouble(dims.z());
                    rb->collisionShapes.push_back(std::make_shared<RBCollisionBox>(rb, p, q, dims));
                }
                break;
            case RBS_COLLISION_SPHERE: {
                    double tx = 0.0, ty = 0.0, tz = 0.0, tr = 0.0;
                    parseDouble(tx);
                    parseDouble(ty);
                    parseDouble(tz);
                    parseDouble(tr);

                    rb->collisionShapes.push_back(std::make_shared<RBCollisionSphere>(rb, P3D(tx, ty, tz), tr));
                }
                break;
            case RBS_COLLISION_PLANE: {
                    double t0 = 0.0, t1 = 0.0, t2 = 0.0, t3 = 0.0, t4 = 0.0, t5 = 0.0;
                    parseDouble(t0); parseDouble(t1); parseDouble(t2);
                    parseDouble(t3); parseDouble(t4); parseDouble(t5);
                    rb->collisionShapes.push_back(std::make_shared<RBCollisionPlane>(rb, P3D(t3, t4, t5), V3D(t0, t1, t2)));
                }
                break;
            case RBS_END_EFFECTOR: {
                    RBEndEffector rbEE;
                    parseDouble(rbEE.pos.x); parseDouble(rbEE.pos.y); parseDouble(rbEE.pos.z);
                    parseDouble(rbEE.radius);
                    rb->endEffectorPoints.push_back(rbEE);
                    rb->collisionShapes.push_back(std::make_shared<RBCollisionSphere>(rb, rbEE.pos, MAX(rbEE.radius, 0.01)));
                } break;
            case RBS_END_RB:
                return rb;  // and... done
                break;
            case RBS_COMMENT:
                parseEntireLine();
                break;
            case RBS_END_OF_FILE:
                break;
            case RBS_UNKNOWN:
                if (verbose) parserError("Unknown token");
                parseEntireLine();
                break;
            default:
                if (verbose) parserError("Unexpected token %s", RBSKeywords[tokenID].keywordText.c_str());
                parseEntireLine();
        }
    }

    throwError("/End_RB expected but not found");
    return rb;
}

RBJoint* RBSLoader::createAndPopulateRBJoint(RBSRepo* rbsRepo) {
    RBJoint* j = new RBJoint();

    while (!endOfTextStream()) {
        int tokenID = getRBSLineType();

        switch (tokenID) {
            case RBS_NAME:
                parseName(j->name);
                break;
            case RBS_PARENT: {
                    std::string tmp;
                    parseName(tmp);
                    j->parent = rbsRepo->getRBByName(tmp.c_str());
                    if (j->parent == nullptr)
                        parserError("Cannot find joint parent RB: %s", tmp.c_str());
                }
                break;
            case RBS_CHILD: {
                    std::string tmp;
                    parseName(tmp);
                    j->child = rbsRepo->getRBByName(tmp.c_str());
                    if (j->child == nullptr)
                        parserError("Cannot find joint child RB: %s", tmp.c_str());
                }
                break;
            case RBS_CPOS:
                parseDouble(j->cJPos.x);
                parseDouble(j->cJPos.y);
                parseDouble(j->cJPos.z);
                break;
            case RBS_PPOS:
                parseDouble(j->pJPos.x);
                parseDouble(j->pJPos.y);
                parseDouble(j->pJPos.z);
                break;

            case RBS_JOINT_AXIS:
                parseDouble(j->rotAxis.x());
                parseDouble(j->rotAxis.y());
                parseDouble(j->rotAxis.z());

                if (j->rotAxis.norm() > 1e-3)
                    j->rotAxis.normalize();
                else
                    j->rotAxis = V3D();
                break;
            case RBS_JOINT_LIMITS:
                parseDouble(j->minAngle);
                parseDouble(j->maxAngle);
                break;
            case RBS_DEFAULT_ANGLE:
                parseDouble(j->defaultJointAngle);
                j->motor.targetPosition = j->defaultJointAngle;
                break;
            case RBS_CONTROL_MODE:
                parseWhitespace();
                if (parseIfKeyword("POS")) j->motor.controlMode = POSITION_CONTROL;
                else if (parseIfKeyword("VEL")) j->motor.controlMode = VELOCITY_CONTROL;
                else if (parseIfKeyword("OFF")) j->motor.controlMode = OFF;
                else if (parseIfKeyword("TORQUE")) j->motor.controlMode = TORQUE_CONTROL;
                else parserError("controlMode: expecting one of POS, VEL, OFF or TORQUE");
                break;
            case RBS_JOINT_MAX_SPEED:
                parseDouble(j->motor.maxSpeed);
                break;
            case RBS_JOINT_MAX_TORQUE:
                parseDouble(j->motor.maxTorque);
                break;
            case RBS_JOINT_END:
                // we now have to link together the child and parent bodies
                if (j->child == nullptr)
                    throwError("Joint %s has no child", j->name.c_str());
                if (j->parent == nullptr)
                    throwError("Joint %s has no parent", j->name.c_str());
                return j;
                break;
            case RBS_COMMENT:
                parseEntireLine();
                break;
            case RBS_UNKNOWN:
                parserError("Unknown token");
                parseEntireLine();
                break;
            case RBS_END_OF_FILE:
                break;
            default:
                parserError("Unexpected token %s", RBSKeywords[tokenID].keywordText.c_str());
                parseEntireLine();
        }
    }

    throwError("/End_Joint expected but not found");
    return j;
}

//i think motor kp/kd and such are written but not read!?

