#ifndef FORCES_H
#define FORCES_H

#include <cmath>
#include <vector>
#include "particle.h"

class Force
{
public:
    Force(void) {}
    virtual ~Force(void) {}

    virtual void apply() = 0;

    void addInfluencedParticle(Particle* p) {
        particles.push_back(p);
    }

    void setInfluencedParticles(const std::vector<Particle*>& vparticles) {
        particles = vparticles;
    }

    void clearInfluencedParticles() {
        particles.clear();
    }

    std::vector<Particle*> getInfluencedParticles() const {
        return particles;
    }

protected:
    std::vector<Particle*>	particles;
    float friction = -0.01f;
};


// Constant acceleration force
class ForceConstAcceleration : public Force
{
public:
    ForceConstAcceleration() { acceleration = Vec3(0,0,0); }
    ForceConstAcceleration(const Vec3& a) { acceleration = a; }
    virtual ~ForceConstAcceleration() {}

    virtual void apply();

    void setAcceleration(const Vec3& a) { acceleration = a; }
    Vec3 getAcceleration() const { return acceleration; }

protected:
    Vec3 acceleration;
};



// Linear drag force
class ForceDragLinear : public Force
{
public:
    ForceDragLinear() { acceleration = Vec3(0,0,0); }
    ForceDragLinear(const Vec3& a) { acceleration = a; }
    virtual ~ForceDragLinear() {}

    virtual void apply();

    void setAcceleration(const Vec3& a) { acceleration = a; }
    Vec3 getAcceleration() const { return acceleration; }

protected:
    Vec3 acceleration;
};



// Quadratic drag force
class ForceDragQuadratic : public Force
{
public:
    ForceDragQuadratic() { acceleration = Vec3(0,0,0); }
    ForceDragQuadratic(const Vec3& a) { acceleration = a; }
    virtual ~ForceDragQuadratic() {}

    virtual void apply();

    void setAcceleration(const Vec3& a) { acceleration = a; }
    Vec3 getAcceleration() const { return acceleration; }

protected:
    Vec3 acceleration;
};


// Black hole gravitational force
class ForceBlackHole: public Force
{
public:
    ForceBlackHole() { position = Vec3(0,0,0); mass = 0.0f;}
    ForceBlackHole(const Vec3& p, const float m) { position = p; mass = m; }
    virtual ~ForceBlackHole() {}

    virtual void apply();

    void setPosition(const Vec3& p) { position = p;}
    Vec3 getPosition() const { return position;}
    void setMass(const float m) { mass = m;}
    float getMass() const { return mass;}

protected:
    Vec3 position;
    float mass;

};


// Spring force
class ForceSpring : public Force
{
public:
    ForceSpring() {}
    ForceSpring(Particle* particle0, Particle* particle1) { p0 = particle0; p1 = particle1; L_distance=(particle1->pos - particle0->pos).norm();}
    virtual ~ForceSpring() {}

    virtual void apply();

public:
    Particle *p0 = nullptr, *p1 = nullptr;

    float ks = 1.f;
    float kd = 0.5f;
    float L_distance;
};




#endif // FORCES_H
