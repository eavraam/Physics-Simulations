#include "model.h"
#include "defines.h"
#include <map>
#include <utility>


Model Model::createQuad() {
    const std::vector<float> verts = {
        -1,0,-1,
        -1,0, 1,
         1,0, 1,
         1,0,-1
    };
    const std::vector<float> norms = {
        0,1,0,
        0,1,0,
        0,1,0,
        0,1,0
    };
    const std::vector<unsigned int> indices = {
        0,1,2,0,2,3
    };

    return Model(verts, norms, indices);
}


Model Model::createCube(bool innerNormals) {
    const std::vector<float> cubeVerts = {
        -1,-1,-1,
        -1,-1, 1,
        -1, 1,-1,
        -1, 1, 1,
         1,-1,-1,
         1,-1, 1,
         1, 1,-1,
         1, 1, 1
    };
    const std::vector<unsigned int> faces = {
         0, 1, 2, 2, 1, 3,
         5, 4, 7, 7, 4, 6,
         0, 4, 1, 1, 4, 5,
         3, 7, 2, 2, 7, 6,
         4, 0, 6, 6, 0, 2,
         1, 5, 3, 3, 5, 7
    };
    const std::vector<float> faceNormals = {
        -1, 0, 0,
         1, 0, 0,
         0,-1, 0,
         0, 1, 0,
         0, 0,-1,
         0, 0, 1
    };

    // repeat vertices to have their normal per face
    std::vector<float> verts(faces.size()*3);
    std::vector<float> norms(faces.size()*3);
    std::vector<unsigned int> indices(faces.size());
    for (unsigned int i = 0; i < faces.size(); i++) {
        unsigned int v = faces[i];
        unsigned int f = i/6;
        verts[3*i  ] = cubeVerts[3*v  ];
        verts[3*i+1] = cubeVerts[3*v+1];
        verts[3*i+2] = cubeVerts[3*v+2];
        norms[3*i  ] = faceNormals[3*f  ] * (innerNormals ? -1 : 1);
        norms[3*i+1] = faceNormals[3*f+1] * (innerNormals ? -1 : 1);
        norms[3*i+2] = faceNormals[3*f+2] * (innerNormals ? -1 : 1);
        indices[i] = i;
    }

    return Model(verts, norms, indices);
}



Model Model::createIcosphere(int numSubdivisions) {

    // adapted from: https://schneide.blog/2016/07/15/generating-an-icosphere-in-c/
    const float X=.525731112119133606f;
    const float Z=.850650808352039932f;
    const float N=0.f;

    std::vector<Vec3> verts = {
        Vec3(-X, N, Z),
        Vec3( X, N, Z),
        Vec3(-X, N,-Z),
        Vec3( X, N,-Z),
        Vec3( N, Z, X),
        Vec3( N, Z,-X),
        Vec3( N,-Z, X),
        Vec3( N,-Z,-X),
        Vec3( Z, X, N),
        Vec3(-Z, X, N),
        Vec3( Z,-X, N),
        Vec3(-Z,-X, N)
    };
    std::vector<unsigned int> tris = {
        0,4,1,  0,9,4,  9,5,4,  4,5,8,  4,8,1,
        8,10,1, 8,3,10, 5,3,8,  5,2,3,  2,7,3,
        7,10,3, 7,6,10, 7,11,6, 11,0,6, 0,1,6,
        6,1,10, 9,0,11, 9,11,2, 9,2,5,  7,2,11
    };

    for (int s = 0; s < numSubdivisions; s++) {

        // keep track of which edges have been already subdivided, since we'll visit them twice
        std::map<std::pair<unsigned int, unsigned int>, unsigned int> edgeVert;

        // for each triangle
        unsigned int currLevelTris = tris.size();
        for (unsigned int t = 0; t < currLevelTris; t += 3) {

            // for each edge between two vertices, create a new vertex if needed
            for (unsigned int i = 0; i < 3; i++) {
                unsigned int v1 = tris[t + i];
                unsigned int v2 = tris[t + (i + 1)%3];
                std::pair<unsigned int, unsigned int> p(v1, v2);

                // create midpoint vertex if not already done
                if (edgeVert.find(p) == edgeVert.end()) {

                    Vec3 p1 = verts[v1];
                    Vec3 p2 = verts[v2];
                    Vec3 newVert = (0.5*(p1 + p2)).normalized();

                    unsigned int idx = verts.size();
                    verts.push_back(newVert);

                    // record vertex index in both senses
                    edgeVert[p] = idx;
                    edgeVert[std::make_pair(v2, v1)] = idx;
                }
            }

            // subdivide triangle into four new triangles (3 new ones, reuse values for central)
            unsigned int v1 = tris[t];
            unsigned int v2 = tris[t+1];
            unsigned int v3 = tris[t+2];
            unsigned int m12 = edgeVert[std::make_pair(v1, v2)];
            unsigned int m23 = edgeVert[std::make_pair(v2, v3)];
            unsigned int m31 = edgeVert[std::make_pair(v3, v1)];
            tris.push_back(v1); tris.push_back(m12); tris.push_back(m31);
            tris.push_back(v2); tris.push_back(m23); tris.push_back(m12);
            tris.push_back(v3); tris.push_back(m31); tris.push_back(m23);
            tris[t] = m12;      tris[t+1] = m23;     tris[t+2] = m31;
        }
    }

    std::vector<float> coords(3*verts.size());
    for (unsigned int i = 0; i < verts.size(); i++) {
        coords[3*i  ] = float(verts[i][0]);
        coords[3*i+1] = float(verts[i][1]);
        coords[3*i+2] = float(verts[i][2]);
    }

    return Model(coords, coords, tris);
}
