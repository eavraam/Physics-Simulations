#include <iostream>
#include "forces.h"

// Constant acceleration force
void ForceConstAcceleration::apply() {
    for (Particle* p : particles) {
        p->force += p->mass * getAcceleration();
    }
}

// Linear drag force
// Fd = -b*v
void ForceDragLinear::apply() {
    for(Particle* p : particles) {
        p->force += friction * p->vel;
    }
}


// Quadratic drag force
// Fd = -c*|v|*v.unit_vector
void ForceDragQuadratic::apply() {
    for(Particle* p : particles) {
        p->force += friction * p->vel * p->vel.norm();
    }
}



// Gravitational force
// Fg = gravity*mass1*mass2 / (distance^2)
void ForceBlackHole::apply() {

    double gravity = 9.81f; // Gravitational constant

    for (Particle* p : particles) {
        double bh_mass = getMass();
        Vec3 bh_pos = getPosition();

        // Calculate the vector from the particle to the black hole
        Vec3 r = bh_pos - p->pos;
        double distance = r.norm();

        // Avoid division by zero and calculate the force
        if (distance > 20.0f) { // 20.0f is a value I found to work "visually nice".
            double force_magnitude = (gravity * bh_mass * p->mass) / (distance * distance);
            Vec3 force = (force_magnitude / (distance*1.5)) * r;
            p->force += force;
        }
    }
}



// Spring force
// if a particle IS fixed, we DO NOT apply forces, else it would be moving.
// if it's NOT fixed, then ADD forces to p0 or SUBTRACT forces from p1.
    // Every particle in a full particle rotation will be once a p0 and once a p1,
        //when computing pairs of particles.
    // If we would add the forces to p0 AND p1, then every particle will ultimately
        // get TWICE the forces.

void ForceSpring::apply() {

    Vec3 deltaPos = p1->pos - p0->pos;
    float distance = deltaPos.norm();
    Vec3 normalizedDeltaPos = deltaPos / distance; // Normalize to get direction

    float F_stretch = ks * (distance - L_distance);
    float F_damp = kd * ((p1->vel - p0->vel).dot(normalizedDeltaPos));

    Vec3 F_total = (F_stretch + F_damp) * normalizedDeltaPos;

//    if (!p0->isFixed()) // I need this on the fixed after all, in case they are released
        p0->force += F_total;
//    if (!p1->isFixed()) // I need this on the fixed after all, in case they are released
        p1->force -= F_total;
}

