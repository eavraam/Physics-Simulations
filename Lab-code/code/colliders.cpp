#include "colliders.h"
#include <cmath>
#include <iostream>


/*
 * Helper functions
 */

void Collider::ApplyElasticCollisionAndFriction(Particle* p, const Vecd& planeN,
                                      double planeD, double kElastic, double kFriction) const
{
    // Elastic collision
    // position: p = p - (1+ε)(np+d)n
    p->pos = p->pos - (1 + kElastic)*(planeN.dot(p->pos) + planeD) * planeN;
    // velocity: v = v - (1+ε)(nv)n = v - (1+ε)v_n
    p->vel = p->vel - (1 + kElastic) * planeN.dot(p->vel) * planeN;

    // Friction forces
    // normal velocity: v_n = (nv)n
    Vecd vel_normal = planeN.dot(p->vel) * planeN;
    // tangent velocity: v_t = v - v_n
    Vecd vel_tangent = p->vel - vel_normal;

    // After applying elasticity and friction
    p->vel = p->vel - kFriction * vel_tangent;
}


/*
 * Plane
 */

bool ColliderPlane::testCollision(Particle* p, double timeStep) const
{
    double current = planeN.dot(p->pos) + planeD;
    double previous = planeN.dot(p->prevPos) + planeD;

    return (current * previous) <= 0.0;
}

void ColliderPlane::resolveCollision(Particle* p, double kElastic, double kFriction) const
{
    // Compute the collision
    ApplyElasticCollisionAndFriction(p, planeN, planeD, kElastic, kFriction);
}

void ColliderPlane::Move(const Vec3& displacement) {}



/*
 * Sphere
 */

bool ColliderSphere::testCollision(Particle* p, double timeStep) const
{
    // (X - C)(X - C)^T = r^2
    // ^T == transposed
    Vecd diff = p->pos - center;
    Vecd diff_transpose = diff.transpose();

    // (diff)(diff_transpose) <= r^2
    if (diff.dot(diff_transpose) <= radius*radius)
    {
        return true;
    }

    return false;
}

void ColliderSphere::resolveCollision(Particle* p, double kElastic, double kFriction) const
{
    // Get the contact point, get the normal from its distance from the center
    // Then: Ax+By+Cz+D=0 => D = - (Ax+By+Cz) (1)
    // We solve so (1) ensures that p->pos lies on the plane defined by planeN, planeD.

    Vecd contactPoint = p->pos - center;
    Vecd planeN = contactPoint.normalized();   

    p->pos = center + planeN * radius;  // fixing the pass-through on large timesteps, potentially introducing some "pop-onto-surface" error.
    double planeD = -planeN.dot(p->pos);

    // Compute the collision
    ApplyElasticCollisionAndFriction(p, planeN, planeD, kElastic, kFriction);
}


void ColliderSphere::Move(const Vec3& displacement)
{
    center += displacement;
}


/*
 * AABB
 */

bool ColliderAABB::testCollision(Particle* p, double timeStep) const
{

    // initialize the aabb intersection limits
    bool xIntersects = ((p->pos.x() + p->radius) >= (aabbPosition.x() - aabbScale.x())) && ((p->pos.x() - p->radius) <= (aabbPosition.x() + aabbScale.x()));
    bool yIntersects = ((p->pos.y() + p->radius) >= (aabbPosition.y() - aabbScale.y())) && ((p->pos.y() - p->radius) <= (aabbPosition.y() + aabbScale.y()));
    bool zIntersects = ((p->pos.z() + p->radius) >= (aabbPosition.z() - aabbScale.z())) && ((p->pos.z() - p->radius) <= (aabbPosition.z() + aabbScale.z()));


    // define the 6 planes of the AABB
    QVector<double> planeD_vector = {
          aabbPosition.x() - aabbScale.x() ,    // xmin
        -(aabbPosition.x() + aabbScale.x()),    // xmax
          aabbPosition.y() - aabbScale.y() ,    // ymin
        -(aabbPosition.y() + aabbScale.y()),    // ymax
          aabbPosition.z() - aabbScale.z() ,    // zmin
        -(aabbPosition.z() + aabbScale.z())     // zmax
    };

    // initialize the needed parameters for the potential collision test
    float CCD_THRESHOLD = timeStep;
    float toi_min = std::numeric_limits<float>::max();
    float vel_component;

    // the particle will eventually hit all 6 planes along its velocity direction
    // so, compute the minimum time of impact (toi), so we know which will be hit first.
    for (unsigned int i = 0; i < 6; i++)
    {
        vel_component = (i % 2 == 0) ? -1.0f : 1.0f; // Extract the current dimension from planeN_vector, in order to handle scalar data.
        float denom = vel_component * p->vel[i / 2]; // Calculate the dot product for this dimension
        if (abs(denom) > 0.0001f)   // some epsilon to avoid div by 0
        {
            float toi = (vel_component * p->pos[i / 2] + planeD_vector[i]) / denom;
            if (toi >= 0)
            {
                if (toi < toi_min)
                {
                    toi_min = toi;
                    p->collision_plane = i;
                }
            }
        }
    }

    // If the minimum toi is less than a continuous collision detection threshold, a collision is possible
    if (toi_min <= CCD_THRESHOLD) {

        // std::cout << "Potential Collision..." << std::endl;

        // If imminent collision ahead (with infinite plane) AND inside the aabb limits
        if (
            (p->collision_plane == 0 && yIntersects && zIntersects ) || // xmin
            (p->collision_plane == 1 && yIntersects && zIntersects ) || // xmax
            (p->collision_plane == 2 && xIntersects && zIntersects ) || // ymin
            (p->collision_plane == 3 && xIntersects && zIntersects ) || // ymax
            (p->collision_plane == 4 && xIntersects && yIntersects ) || // zmin
            (p->collision_plane == 5 && xIntersects && yIntersects )    // zmax
            )
        {
            return true; // There is collision with the aabb.
        }

    }

    return false;   // No collision.
}


void ColliderAABB::resolveCollision(Particle* p, double kElastic, double kFriction) const
{
    // std::cout << "---------- Collided ----------" << std::endl;

    // define the planeD for each of the 6 axes
    double planeD;
    switch (p->collision_plane)
    {
        case 0: // xmin
            planeD = aabbPosition.x() - aabbScale.x();
            break;
        case 1: // xmax
            planeD = -(aabbPosition.x() + aabbScale.x());
            break;
        case 2: // ymin
            planeD = aabbPosition.y() - aabbScale.y();
            break;
        case 3: // ymax
            planeD = -(aabbPosition.y() + aabbScale.y());
            break;
        case 4: // zmin
            planeD = aabbPosition.z() - aabbScale.z();
            break;
        case 5: // zmax
            planeD = -(aabbPosition.z() + aabbScale.z());
            break;
    }

    // Compute the collision
    ApplyElasticCollisionAndFriction(p, planeN_vector[p->collision_plane], planeD, kElastic, kFriction);

}

void ColliderAABB::Move(const Vec3& displacement)
{
    aabbPosition += displacement;
}


/*
* Particles
*/
// Collisions between Particles
// Similar to ColliderSphere, modified for particles
bool testParticleCollision(Particle* p1, Particle* p2)
{
    Vecd diff = p1->pos - p2->pos;
    Vecd diff_transpose = diff.transpose();

    // (diff)(diff_transpose) <= r^2
    // particles have same radius
    if (diff.dot(diff_transpose) <= p1->radius*p2->radius)
    {
        return true;
    }

    return false;
}

// Resolve particle-particle collision
void resolveParticleCollision(Particle* p1, Particle* p2, double kElastic, double kFriction)
{
    // Get the contact point, get the normal from its distance from the center
    // Then: Ax+By+Cz+D=0 => D = - (Ax+By+Cz) (1)
    // We solve so (1) ensures that p->pos lies on the plane defined by planeN, planeD.
    Vecd contactPoint = p1->pos - p2->pos;
    Vecd planeN = contactPoint.normalized();

    p1->pos = p2->pos + planeN * p2->radius;  // fixing the pass-through on large timesteps, potentially introducing some "pop-onto-surface" error.
    double planeD = -planeN.dot(p1->pos);

    // Elastic collision
    // position: p = p - (1+ε)(np+d)n
    p1->pos = p1->pos - (1 + kElastic)*(planeN.dot(p1->pos) + planeD) * planeN;
    // velocity: v = v - (1+ε)(nv)n = v - (1+ε)v_n
    p1->vel = p1->vel - (1 + kElastic) * planeN.dot(p1->vel) * planeN;

    // Friction forces
    // normal velocity: v_n = (nv)n
    Vecd vel_normal = planeN.dot(p1->vel) * planeN;
    // tangent velocity: v_t = v - v_n
    Vecd vel_tangent = p1->vel - vel_normal;

    // After applying elasticity and friction
    p1->vel = p1->vel - kFriction * vel_tangent;
}
