#pragma once

#include <gui/guiMath.h>
#include <gui/shader.h>

// possible loss of data in conversion between double and float
#pragma warning(disable : 4244)
// deprecated/unsafe functions such as fopen
#pragma warning(disable : 4996)

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/hash.hpp>
#include <map>
#include <vector>


struct Vertex {
    glm::vec3 position = glm::vec3(0, 0, 0);
    glm::vec3 normal = glm::vec3(0, 0, 0);
    glm::vec2 texCoords = glm::vec2(0, 0);

    Vertex() {}
    Vertex(glm::vec3 position, glm::vec3 normal, glm::vec2 texCoords)
        : position(position), normal(normal), texCoords(texCoords) {}

    bool operator==(const Vertex &other) const {
        return position == other.position && normal == other.normal &&
               texCoords == other.texCoords;
    }
};

struct Texture {
    unsigned int id;
    std::string path;
};

class Mesh {
public:
    enum TextureType { DIFFUSE, SPECULAR, NORMAL, AMBIENT };
    typedef std::map<TextureType, std::vector<Texture>> TextureMap;

    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;
    TextureMap textures;

    bool showTextures = true;

private:
    unsigned int VAO, VBO, EBO;

public:
    Mesh(std::vector<Vertex> vertices, std::vector<unsigned int> indices,
         std::map<TextureType, std::vector<Texture>> textures);

    Mesh(std::vector<Vertex> vertices, std::vector<unsigned int> indices);

    /**
     * render the mesh using shader and specified color
     */
    void draw(const V3D &color, double alpha = 1.0) const;

    /**
     * render the mesh using shader and default color
     */
    void draw(double alpha = 1.0) const;

    void reinitialize(std::vector<Vertex> vertices,
                      std::vector<unsigned int> indices);

private:
    // initializes all the buffer objects/arrays
    void setupMesh();
};

