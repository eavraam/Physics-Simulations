#ifndef SCENEFLUID_H
#define SCENEFLUID_H

#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <list>
#include "scene.h"
#include "widgetfluid.h"
#include "particlesystem.h"
#include "integrators.h"
#include "colliders.h"
#include "spatialhash.h"

class SceneFluid : public Scene
{
    Q_OBJECT

public:
    SceneFluid();
    virtual ~SceneFluid();

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

    // Check intersections
    bool checkIntersectionWithAABB(const Vec3& rayOrigin, const Vec3& rayDir, const ColliderAABB& collider);
    bool checkIntersectionWithSphere(const Vec3& rayOrigin, const Vec3& rayDir, const ColliderSphere& collider);

    // Create particles
    void createBoundaryParticles();
    void createWaterParticles();

    // SPH-related methods
    void updateSPH();
    double calculateDensity(Particle* p);
    double calculatePressure(double density);
    Vec3 calculateForces(Particle* p);



public slots:
    void updateSimParams();

protected:
    WidgetFluid* widget = nullptr;

    QOpenGLShaderProgram* shader            = nullptr;
    QOpenGLVertexArrayObject* vaoSphereS    = nullptr;
    QOpenGLVertexArrayObject* vaoFloor      = nullptr;
    QOpenGLVertexArrayObject* vaoWall      = nullptr;
    QOpenGLVertexArrayObject* vaoPropSphere = nullptr;
    QOpenGLVertexArrayObject* vaoPropBox    = nullptr;
    QOpenGLVertexArrayObject* vaoBlackHole  = nullptr;
    unsigned int numFacesSphereS    = 0;
    unsigned int numFacesPropSphere = 0;
    unsigned int numFacesPropBox    = 0;
    unsigned int numFacesBlackHole  = 0;

    IntegratorSymplecticEuler integrator;
    ParticleSystem system;
    std::list<Particle*> deadParticles;

    ForceConstAcceleration* fGravity;
    ForceDragLinear* fDragLinear;
    ForceDragQuadratic* fDragQuadratic;
    ForceBlackHole* fBlackHole;

    ColliderPlane colliderFloor;
    ColliderAABB colliderPropBox;
    std::vector<Collider*> colliders;
    Collider* collidedCollider = nullptr;
    // testing for better boundaries
    ColliderPlane colliderCeiling;
    ColliderPlane colliderLeftWall;
    ColliderPlane colliderRightWall;
    ColliderPlane colliderFrontWall;
    ColliderPlane colliderBackWall;
    // end of testing

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

    QVector<Particle*> boundaryParticles;
    QVector<Particle*> waterParticles;

    Vec3 boundaryArea = Vec3(40, 80, 40);

    // SPH parameters
    double restDensity = 1000.0f;
    double k_const = 1.0f;
    double viscosityConstant = 0.01f;
    double kernelRadius; // the smoothing radius for SPH kernel

    int mouseX, mouseY;
};

#endif // SCENEFLUID_H

