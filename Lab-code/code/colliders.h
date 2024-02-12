#ifndef COLLIDERS_H
#define COLLIDERS_H

#include <QVector>
#include <limits>
#include "defines.h"
#include "particle.h"


// Abstract interface
class Collider
{
public:
    Collider() {}
    virtual ~Collider() {}

    virtual void ApplyElasticCollisionAndFriction(Particle* p, const Vecd& planeN, double planeD, double kElastic, double kFriction) const;

    virtual bool testCollision(Particle* p, double timeStep) const = 0;
    virtual void resolveCollision(Particle* p, double kElastic, double kFriction) const = 0;
    virtual void Move(const Vec3& displacement) = 0;
};



// Plane collider
class ColliderPlane : public Collider
{
public:
    ColliderPlane() { planeN = Vec3(0,0,0); planeD = 0; }
    ColliderPlane(const Vec3& n, double d) : planeN(n), planeD(d) {}
    virtual ~ColliderPlane() {}

    void setPlane(const Vec3& n, double d) { this->planeN = n; this->planeD = d; }

    virtual bool testCollision(Particle* p, double timeStep) const;
    virtual void resolveCollision(Particle* p, double kElastic, double kFriction) const;
    virtual void Move(const Vec3& displacement);

protected:
    Vec3 planeN;
    double planeD;
};



// Sphere collider
class ColliderSphere : public Collider
{
public:
    ColliderSphere() { center = Vec3(0,0,0); radius = 0; }
    ColliderSphere(const Vec3& c, const double r) : center(c), radius(r) {}
    virtual ~ColliderSphere() {}

    void setSphere(const Vec3& c, const double r) { this->center = c; this->radius = r; }

    virtual bool testCollision(Particle* p,double timeStep) const;
    virtual void resolveCollision(Particle* p, double kElastic, double kFriction) const;
    virtual void Move(const Vec3& displacement);

    Vec3 getCenter() const { return center; }
    double getRadius() const { return radius; }

public :
    Vec3 center;
    double radius;
};

// AABB collider
class ColliderAABB : public Collider
{
public:
    ColliderAABB() { aabbPosition = Vec3(0,0,0); aabbScale = Vec3(0,0,0); }
    ColliderAABB(const Vec3& p, const Vec3& s) : aabbPosition(p), aabbScale(s) {}
    virtual ~ColliderAABB() {}

    void setAABB(const Vec3& p, const Vec3& s) { this->aabbPosition = p; this->aabbScale = s; }

    virtual bool testCollision(Particle* p, double timeStep) const;
    virtual void resolveCollision(Particle* p, double kElastic, double kFriction) const;
    virtual void Move(const Vec3& displacement);

public:
    Vec3 aabbPosition, aabbScale;

protected:
    QVector<Vec3> planeN_vector = {
        Vec3(-1.0f,  0.0f,  0.0f),  // xmin
        Vec3( 1.0f,  0.0f,  0.0f),  // xmax
        Vec3( 0.0f, -1.0f,  0.0f),  // ymin
        Vec3( 0.0f,  1.0f,  0.0f),  // ymax
        Vec3( 0.0f,  0.0f, -1.0f),  // zmin
        Vec3( 0.0f,  0.0f,  1.0f)   // zmax
    };
};


// Particle-particle collisions
bool testParticleCollision(Particle* p1, Particle* p2);
void resolveParticleCollision(Particle* p1, Particle* p2, double kElastic, double kFriction);


#endif // COLLIDERS_H
