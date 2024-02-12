#ifndef MODEL_H
#define MODEL_H

#include <vector>

/*
 *  This class represents an indexed triangle mesh model with normals per vertex.
 *  If you want normals per face, vertices need to be repeated for each face with the corresponding normal.
 */
class Model
{
public:

    Model() {}
    Model(const std::vector<float>& v, const std::vector<float>& n, const std::vector<unsigned int>& i) :
        vertexCoords(v), vertexNormals(n), vertexIndices(i) {}
    virtual ~Model() {}

    int numVertices() {
        return vertexCoords.size()/3;
    }
    int numFaces() {
        return vertexIndices.size()/vertsPerFace;
    }
    int verticesPerFace() {
        return vertsPerFace;
    }

    const std::vector<float>& getVertexCoords() { return vertexCoords; }
    const std::vector<float>& getNormals() { return vertexNormals; }
    const std::vector<unsigned int>& getIndices() { return vertexIndices; }

    const float* getCoordsPtr() { return &vertexCoords[0]; }
    const float* getNormalsPtr() { return &vertexNormals[0]; }
    const unsigned int* getIndicesPtr() { return &vertexIndices[0]; }

protected:

    std::vector<float> vertexCoords;
    std::vector<float> vertexNormals;
    std::vector<unsigned int> vertexIndices;
    static constexpr int vertsPerFace = 3;


// some static Model constructors for convenience
public:
    // quad from (-1,-1,0) to (1,1,0)
    static Model createQuad();
    // cube form (-1,-1,-1) to (1,1,1)
    static Model createCube(bool innerNormals=false);
    // radius 1 icosphere created from numSubdivisions of an icosahedron. 0 returns icosahedron.
    static Model createIcosphere(int numSubdivisions = 0);

};



#endif // MODEL_H
