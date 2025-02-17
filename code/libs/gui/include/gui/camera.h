#pragma once

#include <glad/glad.h>
#include <gui/guiMath.h>
#include <gui/inputstate.h>

#include <utils/geoms.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>
#include <vector>


//implements a simple tracking camera
class Camera {
public:
    // Camera Attributes
    // explicitely keep track of the rotations about the "up" and "horizontal
    // right" directions to avoid singularities in decomposition of large
    // rotations
    float rotAboutUpAxis = 0, rotAboutRightAxis = 0;
    // the location of the camera is determined as a function of its
    // orientation, its target, and how far along the target it is
    float distanceToTarget = 50;
    // this is the target of the camera, in world coordinates
    glm::vec3 target = {0, 0, 0};
    // we need to keep track of the view direction (default down the z axis) and
    // up axis
    glm::vec3 direction = {0, 0, -1}, up = {0, 1, 0};

    bool freezeCameraOrientation = false;

    // params to keep track of perspective projection map
    float aspectRatio = 1.0;
    float fieldOfView = 45.0;
    float zNear = 0.1f, zFar = 1000.0f;

    // params that control the camera movement
    float mouseSensitivity = 0.01f;

    // Constructor with vectors
    Camera(float distToTarget = 50) { distanceToTarget = distToTarget; }

    // Returns the view matrix calculated using Euler Angles and the LookAt
    // Matrix
    glm::mat4 getViewMatrix() const {
        glm::vec3 worldUp = rotation() * glm::vec4(up, 1);

        return glm::lookAt(position(), target, worldUp);
    }

    glm::vec3 position() const {
        return target + glm::vec3(rotation() * glm::vec4(direction, 1)) *
                            -distanceToTarget;
    }

    // Processes input received from a mouse input system. Expects the offset
    // value in both the x and y direction.
    void processLeftMouseMovement(float xOffset, float yOffset) {
        if (freezeCameraOrientation) return;
        rotAboutUpAxis += xOffset * mouseSensitivity;
        rotAboutRightAxis += yOffset * mouseSensitivity;
        if (rotAboutRightAxis > 1.5) rotAboutRightAxis = 1.5;
        if (rotAboutRightAxis < -1.5) rotAboutRightAxis = -1.5;
    }

    void processRightMouseMovement(float xOffset, float yOffset) {
        glm::vec3 m = {xOffset, yOffset, 0.f};
        glm::mat4 rot = getViewMatrix();
        target += 0.0005f * glm::vec3(glm::inverse(rot) * glm::vec4(m, 0.f)) *
                  distanceToTarget;
    }

    // Processes input received from a mouse scroll-wheel event. Only requires
    // input on the vertical wheel-axis
    void processMouseScroll(float xOffset, float yOffset) {
        distanceToTarget *= 1.0f - yOffset * 0.05f;
        if (distanceToTarget < 0.1) distanceToTarget = 0.1f;
    }

    // Processes input received from a mouse scroll-wheel event. Only requires
    // input on the vertical wheel-axis
    void processMouseMove(MouseState ms, KeyboardState ks) {
        if (ms.dragging == false) return;
        if (ms.lButtonPressed && ks[GLFW_KEY_LEFT_ALT] == false)
            processLeftMouseMovement((float)ms.mouseMoveX,
                                     (float)ms.mouseMoveY);
        if (ms.rButtonPressed ||
            (ms.lButtonPressed && ks[GLFW_KEY_LEFT_ALT] == true))
            processRightMouseMovement((float)ms.mouseMoveX,
                                      (float)ms.mouseMoveY);
    }

    glm::mat4 getProjectionMatrix() const {
        return glm::perspective(glm::radians(fieldOfView), aspectRatio, zNear,
                                zFar);
    }

    Ray getRayFromScreenCoordinates(double xpos, double ypos){
        P3D rayPos;
        V3D rayDir;
        getRayFromScreenCoordinates(xpos, ypos, rayPos, rayDir);
        return Ray(rayPos, rayDir);
    }

    // xpos and ypos represent points in screen coordinates (e.g. where the
    // mouse cursor is); rayOrigin and rayDirection will represent quantities in
    // world coordinates...
    void getRayFromScreenCoordinates(double xpos, double ypos, P3D &rayOrigin,
                                     V3D &rayDirection) {
        // this method assumes the correct viewport is set when this method is
        // called - can get slightly tricky when multiple sub windows are
        // used...
        glm::vec4 viewportParams;
        glGetFloatv(GL_VIEWPORT, &viewportParams[0]);

        // 0 and 1 are here the near and far z planes in normalized coordinates
        // (e.g. right before projection to screen coordinates)
        glm::vec3 p1 = glm::unProject(
            glm::vec3(xpos, viewportParams[3] - ypos - 1, 0), getViewMatrix(),
            getProjectionMatrix(), viewportParams);
        glm::vec3 p2 = glm::unProject(
            glm::vec3(xpos, viewportParams[3] - ypos - 1, 1), getViewMatrix(),
            getProjectionMatrix(), viewportParams);

        rayOrigin = toP3D(p1);
        rayDirection = toV3D(p2 - p1);
        rayDirection.normalize();
    }

    glm::mat4 rotation() const {
        glm::mat4 rot(1.f);
        rot = glm::rotate(rot, rotAboutUpAxis, up);
        rot = glm::rotate(rot, rotAboutRightAxis, glm::cross(up, direction));
        return rot;
    }

    void setRotation(const glm::mat4& rot){
        setRotation(glm::quat_cast(rot));
    }

    void setRotation(const glm::quat& q){
        glm::vec3 eulerAngles = glm::eulerAngles(q); //euler angles yxz

        Quaternion qQ = toQuaternion(q);
        double tmp1, tmp2, tmp3;
        qQ.computeEulerAngles(toV3D(direction), toV3D(glm::cross(up, direction)), toV3D(up), tmp1, tmp2, tmp3);

        clamp(&tmp2, -1.5, 1.5);

        rotAboutRightAxis = tmp2;
        rotAboutUpAxis = tmp3;
    }

};

