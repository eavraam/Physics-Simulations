#ifndef SCENECLOTH_H
#define SCENECLOTH_H

#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include "scene.h"
#include "widgetcloth.h"
#include "particlesystem.h"
#include "integrators.h"
#include "colliders.h"
#include "spatialhash.h"


class SceneCloth : public Scene
{
    Q_OBJECT

public:
    SceneCloth();
    virtual ~SceneCloth();

    virtual void initialize(double dt, int dr, double kEl, double kFr, int bht);
    virtual void reset(double dt, int dr, double kEl, double kFr, int bht);
    virtual void update(double dt);
    virtual void paint(const Camera& cam);

    virtual void mousePressed(const QMouseEvent* e, const Camera& cam);
    virtual void mouseMoved(const QMouseEvent* e, const Camera& cam);
    virtual void mouseReleased(const QMouseEvent* e, const Camera& cam);
    virtual void keyPressed(const QKeyEvent* e, const Camera& cam);

    virtual void getSceneBounds(Vec3& bmin, Vec3& bmax) {
        bmin = Vec3(-100, -100, -100);
        bmax = Vec3( 100,  100,  100);
    }
    virtual unsigned int getNumParticles() { return system.getNumParticles(); }

    virtual QWidget* sceneUI() { return widget; }

    bool checkIntersectionWithParticle(const Vec3& rayOrigin, const Vec3& rayDir, Particle* p);
    bool checkIntersectionWithAABB(const Vec3& rayOrigin, const Vec3& rayDir, const ColliderAABB& collider);
    bool checkIntersectionWithSphere(const Vec3& rayOrigin, const Vec3& rayDir, const ColliderSphere& collider);

    void computeCloth(int particles_width, int particles_height,
                            QVector<Particle*>& sceneParticles, QVector<ForceSpring*>& sceneForceSprings,
                            double clothWidth, double clothHeight, double particleRadius);

public slots:
    void updateSprings();
    void updateSimParams();
    void freeAnchors();

protected:
    // ui
    WidgetCloth* widget = nullptr;

    // opengl & render
    QOpenGLShaderProgram* shaderPhong = nullptr;
    QOpenGLShaderProgram* shaderCloth = nullptr;
    QOpenGLVertexArrayObject* vaoSphereS = nullptr;
    QOpenGLVertexArrayObject* vaoSphereL = nullptr;
    QOpenGLVertexArrayObject* vaoCube    = nullptr;
    QOpenGLVertexArrayObject* vaoMesh    = nullptr;
    QOpenGLVertexArrayObject* vaoFloor      = nullptr;
    QOpenGLVertexArrayObject* vaoPropSphere = nullptr;
    QOpenGLVertexArrayObject* vaoPropBox    = nullptr;
    QOpenGLVertexArrayObject* vaoBlackHole    = nullptr;
    QOpenGLBuffer* vboMesh = nullptr;
    QOpenGLBuffer* iboMesh = nullptr;
    unsigned int numFacesSphereS = 0, numFacesSphereL = 0;
    unsigned int numMeshIndices = 0;
    unsigned int numFacesPropSphere = 0;
    unsigned int numFacesPropBox    = 0;
    unsigned int numFacesBlackHole  = 0;
    bool showParticles = true;

    double timeStep;
    int drag;
    double kElastic, kFriction;
    double emitRate;

    // Colliders
    ColliderPlane colliderFloor;
    ColliderSphere colliderPropSphere;
    ColliderAABB colliderPropBox;

    // physics
    IntegratorVerlet integrator; // Verlet is the way to go
    ParticleSystem system;
    ForceConstAcceleration* fGravity = nullptr;
    ForceDragLinear* fDragLinear;
    ForceDragQuadratic* fDragQuadratic;
    ForceBlackHole* fBlackHole = nullptr;

    QVector<Particle*> sceneParticles;
    QVector<ForceSpring*> sceneForceSprings;
    SpatialHash* spatialHash;

    int blackHoleToggle;
    Vec3 blackHolePos;
    float blackHoleRadius;
    float blackHoleMass;

    // cloth properties
    std::vector<bool> fixedParticle;
    double clothWidth, clothHeight;
    int numParticles, numParticlesX=60, numParticlesY=100;
    int selectedParticle = -1;
    const int relaxationSteps = 10;
    bool spatial_hash_toggle = false;

    // collision properties
    bool checkCollisions = true;
    double colBounce = 0.01;
    double colFriction = 0.05;
    double particleRadius = 1;
    float collisionRadius = 2.5f; // Slightly more than twice the particle radius

    // mouse interaction
    int grabX, grabY;
    Vec3 cursorWorldPos;

//    bool testParticleCollision(Particle* p1, Particle* p2) {
//        Vec3 diff = p1->pos - p2->pos;
//        double distanceSquared = diff.squaredNorm();
//        double radiusSum = p1->radius + p2->radius;
//        return distanceSquared <= (radiusSum * radiusSum);
//    }

//    void resolveParticleCollision(Particle* p1, Particle* p2) {
//        Vec3 midPoint = (p1->pos + p2->pos) / 2.0;
//        Vec3 direction = (p1->pos - p2->pos).normalized();

//        double overlap = 2.0 * p1->radius - (midPoint - p1->pos).norm();
//        p1->pos = midPoint + direction * overlap / 2.0;
//        p2->pos = midPoint - direction * overlap / 2.0;
//    }



};

#endif // SCENECLOTH_H
