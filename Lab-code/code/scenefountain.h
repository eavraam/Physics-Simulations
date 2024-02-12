#ifndef SCENEFOUNTAIN_H
#define SCENEFOUNTAIN_H

#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <list>
#include "scene.h"
#include "widgetfountain.h"
#include "particlesystem.h"
#include "integrators.h"
#include "colliders.h"
#include "spatialhash.h"

class SceneFountain : public Scene
{
    Q_OBJECT

public:
    SceneFountain();
    virtual ~SceneFountain();

    //virtual void initialize(int dr, double kEl, double kFr);
    virtual void initialize(double dt, int dr, double kEl, double kFr, int bht);
    //virtual void reset(int dr, double kEl, double kFr);
    virtual void reset(double dt, int dr, double kEl, double kFr, int bht);
    virtual void update(double dt);
    virtual void paint(const Camera& cam);

    virtual void mousePressed(const QMouseEvent* e, const Camera& cam);
    virtual void mouseMoved(const QMouseEvent* e, const Camera& cam);

    virtual void getSceneBounds(Vec3& bmin, Vec3& bmax) {
        bmin = Vec3(-50, -10, -50);
        bmax = Vec3( 50, 100, 50);
    }
    virtual unsigned int getNumParticles() { return system.getNumParticles(); }

    virtual QWidget* sceneUI() { return widget; }

    bool checkIntersectionWithAABB(const Vec3& rayOrigin, const Vec3& rayDir, const ColliderAABB& collider);
    bool checkIntersectionWithSphere(const Vec3& rayOrigin, const Vec3& rayDir, const ColliderSphere& collider);

public slots:
    void updateSimParams();

protected:
    WidgetFountain* widget = nullptr;

    QOpenGLShaderProgram* shader            = nullptr;
    QOpenGLVertexArrayObject* vaoSphereS    = nullptr;
    QOpenGLVertexArrayObject* vaoFloor      = nullptr;
    QOpenGLVertexArrayObject* vaoPropSphere = nullptr;
    QOpenGLVertexArrayObject* vaoPropBox    = nullptr;
    QOpenGLVertexArrayObject* vaoBlackHole  = nullptr;
    unsigned int numFacesSphereS    = 0;
    unsigned int numFacesPropSphere = 0;
    unsigned int numFacesPropBox    = 0;
    unsigned int numFacesBlackHole  = 0;

    IntegratorRK2 integrator;
    ParticleSystem system;
    std::list<Particle*> deadParticles;

    ForceConstAcceleration* fGravity;
    ForceDragLinear* fDragLinear;
    ForceDragQuadratic* fDragQuadratic;
    ForceBlackHole* fBlackHole;

    ColliderPlane colliderFloor;
    ColliderSphere colliderPropSphere;
    ColliderAABB colliderPropBox;
    std::vector<Collider*> colliders;
    Collider* collidedCollider = nullptr;

    double timeStep;
    int drag;
    double kElastic, kFriction;
    double emitRate;
    double maxParticleLife;

    SpatialHash *spatialHash = nullptr;

    int blackHoleToggle;
    Vec3 blackHolePos;
    float blackHoleRadius;
    float blackHoleMass;

    Vec3 fountainPos;
    int mouseX, mouseY;
};

#endif // SCENEFOUNTAIN_H
