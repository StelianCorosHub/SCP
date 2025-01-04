#pragma once

#include <gui/guiMath.h>
#include <gui/mesh.h>
#include <gui/shader.h>

//TODO: we need the shaders to have a callback function where they set their own properties according to the
// model/mesh we're trying to draw. At the moment, the implementation makes too many assupmtions
// (e.g. all shaders have support for textures).

unsigned int textureFromFile(const char *path, const std::string &directory);

class Model {
public:
    std::vector<Mesh::TextureMap>
        material_textures;  // stores all the textures loaded so far,
                            // optimization to make sure textures aren't loaded
                            // more than once.
    std::vector<Mesh> meshes;

    // this is the name of the model, in case it was loaded from a file
    std::string mName;

    // scale about the x, y, and z axes, applied before any other
    // transformations
    V3D scale = V3D(1, 1, 1);
    // the position of the model in world coordinates
    P3D position = P3D(0, 0, 0);
    // the orientation of the model - takes vectors from local coordinates to
    // world coordinates
    Quaternion orientation = Quaternion::Identity();

    bool highlighted = false;
    bool selected = false;

    Model();
    Model(std::string const &path);

    void draw(const V3D &color,
              const glm::mat4 &transform, double alpha = 1.0) const;
    void draw(const V3D &color, double alpha = 1.0) const;
    void draw(double alpha = 1.0) const;

    void setDisplayTexturesFlag(bool showTextures);

    glm::mat4 getTransform() const;

private:
    // loads a model with tinyobjloader from file and stores the resulting
    // meshes in the meshes vector.
    void loadModel(std::string const &path_);

public:
    bool hitByRay(const P3D &r_o, const V3D &r_v, P3D &hitPoint, double &t,
                  V3D &n) const;

    bool hitByRay(const P3D &r_o, const V3D &r_v) const;

    bool hitByRay(const P3D &r_o, const V3D &r_v, double &t) const;

    bool hitByRay(const P3D &r_o, const V3D &r_v, P3D &hitPoint) const;

    bool hitByRay(const P3D &r_o, const V3D &r_v, P3D &hitPoint,
                  V3D &hitNormal) const;
};

inline unsigned int textureFromFile(const char *path,
                                    const std::string &directory);

inline Model getGroundModel(double s = 100) {
    std::vector<Vertex> vertices = {
            {glm::vec3(-s, 0, -s), glm::vec3(0, 1, 0), glm::vec2(0, 0)},
            {glm::vec3(s, 0, -s), glm::vec3(0, 1, 0), glm::vec2(1, 0)},
            {glm::vec3(s, 0, s), glm::vec3(0, 1, 0), glm::vec2(1, 1)},
            {glm::vec3(-s, 0, s), glm::vec3(0, 1, 0), glm::vec2(0, 1)}
    };

    std::vector<unsigned int> indices = {0, 2, 1, 0, 3, 2};
    Mesh groundMesh(vertices, indices);
    Model ground;
    ground.meshes.push_back(groundMesh);

    return ground;
}

class SimpleGroundModel {
public:
    Model ground = getGroundModel(20);
    Model grid1 = Model(SCP_DATA_FOLDER "/meshes/grid1.obj");
    Model grid2 = Model(SCP_DATA_FOLDER "/meshes/grid2.obj");

    void draw(const V3D &col = V3D(0.7, 0.7, 0.9)) {
        grid1.draw(V3D(0.1, 0.1, 0.1));
        grid2.draw(V3D(0.5, 0.5, 0.5));
        ground.draw(col);
    }
};

