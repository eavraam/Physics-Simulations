#include "scenefluid.h"
#include "glutils.h"
#include "model.h"
#include <QOpenGLFunctions_3_3_Core>


SceneFluid::SceneFluid() {
    widget = new WidgetFluid();
    connect(widget, SIGNAL(updatedParameters()), this, SLOT(updateSimParams()));
}


SceneFluid::~SceneFluid() {
    if (widget)         delete widget;
    if (shader)         delete shader;
    if (fGravity)       delete fGravity;
    if (fBlackHole)     delete fBlackHole;
    if (vaoFloor)       delete vaoFloor;
    if (vaoSphereS)     delete vaoSphereS;
    if (vaoPropSphere)  delete vaoPropSphere;
    if (vaoPropBox)     delete vaoPropBox;
    if (vaoBlackHole)   delete vaoBlackHole;
}


void SceneFluid::initialize(double dt, int dr, double kEl, double kFr, int bht)
{
    // Initialize drag, elasticity, friction and black hole
    timeStep = dt;
    drag = dr;
    kElastic = kEl;
    kFriction = kFr;
    blackHoleToggle = bht;

    // load shader
    shader = glutils::loadShaderProgram(":/shaders/phong.vert", ":/shaders/phong.frag");

    // create floor VAO
    Model quad = Model::createQuad();
    vaoFloor = glutils::createVAO(shader, &quad);
    glutils::checkGLError();

    // create sphere VAO
    Model sphere = Model::createIcosphere(1);
    vaoSphereS = glutils::createVAO(shader, &sphere);
    numFacesSphereS = sphere.numFaces();
    glutils::checkGLError();

    // create props' VAOs
    Model propSphere = Model::createIcosphere(1);
    vaoPropSphere = glutils::createVAO(shader, &propSphere);
    numFacesPropSphere = propSphere.numFaces();
    glutils::checkGLError();

    Model propBox = Model::createCube();
    vaoPropBox = glutils::createVAO(shader, &propBox);
    numFacesPropBox = propBox.numFaces();
    glutils::checkGLError();

    // create black hole VAO
    Model blackHole = Model::createIcosphere(1);
    vaoBlackHole = glutils::createVAO(shader, &blackHole);
    numFacesBlackHole = blackHole.numFaces();
    glutils::checkGLError();


    // Forces
    // we need to create one per system to assign its particle
    // reset them in the init, just in case
    system.clearForces();

    // gravity
    fGravity = new ForceConstAcceleration();
    system.addForce(fGravity); // always applied

    // drag (linear, quadratic)
    fDragLinear = new ForceDragLinear();
    fDragQuadratic = new ForceDragQuadratic();
    if (drag == 1)
    {
        system.addForce(fDragLinear);
    }
    else if (drag == 2)
    {
        system.addForce(fDragQuadratic);
    }

    // black hole
    blackHolePos = Vec3(-80.0f, 100.0f, -10.0f);
    blackHoleRadius = 5.0f;
    blackHoleMass = 10000.0f;

    fBlackHole = new ForceBlackHole(blackHolePos, blackHoleMass);
    if (blackHoleToggle)
    {
        system.addForce(fBlackHole);
    }
}


void SceneFluid::reset(double dt, int dr, double kEl, double kFr, int bht)
{

    // Reset drag, elasticity, friction and black hole
    timeStep = dt;
    drag = dr;
    kElastic = kEl;
    kFriction = kFr;
    blackHoleToggle = bht;

    std::cout << "Reset kElastic to: " << kElastic << std::endl;
    std::cout << "Reset kFriction to: " << kFriction << std::endl;
    std::cout << "Reset drag to: " << drag << std::endl;
    std::cout << "Reset black hole toggle to: " << blackHoleToggle << std::endl;

    // reset random seed
    Random::seed(1337);

    // scene
    colliderPropBox.setAABB(Vec3(10.0f, 60.0f, 10.0f), boundaryArea);

    // Floor and Ceiling
    colliderFloor.setPlane(Vec3(0, 1, 0), 0);
    colliderCeiling.setPlane(Vec3(0, -1, 0), boundaryArea.y());

    // Walls
    colliderLeftWall.setPlane(Vec3(1, 0, 0), boundaryArea.x());
    colliderRightWall.setPlane(Vec3(-1, 0, 0), boundaryArea.x());
    colliderFrontWall.setPlane(Vec3(0, 0, -1), boundaryArea.z());
    colliderBackWall.setPlane(Vec3(0, 0, 1), boundaryArea.z());


    // Reset particles and forces
    system.clearForces();
    system.deleteParticles();
    boundaryParticles.clear();
    waterParticles.clear();

    // clear all particle forces
    // NULL pointer reference, need to handle it.
    if (fGravity != nullptr)
        fGravity->clearInfluencedParticles();
    if (fDragLinear != nullptr)
        fDragLinear->clearInfluencedParticles();
    if (fDragQuadratic != nullptr)
        fDragQuadratic->clearInfluencedParticles();
    if (fBlackHole != nullptr)
        fBlackHole->clearInfluencedParticles();

    // Forces
    fGravity = new ForceConstAcceleration();
    system.addForce(fGravity);

    if (blackHoleToggle)
    {
        fBlackHole = new ForceBlackHole(blackHolePos, blackHoleMass);
        system.addForce(fBlackHole);
    }
    if (drag == 1)
    {
        fDragLinear = new ForceDragLinear();
        system.addForce(fDragLinear);
    }
    else if (drag == 2)
    {
        fDragQuadratic = new ForceDragQuadratic();
        system.addForce(fDragQuadratic);
    }

    // Create boundary particles
//    createBoundaryParticles();

    // Create water particles
    createWaterParticles();

    // Add gravity and black hole forces to the particles
    for (Particle* p : waterParticles)
    {
        fGravity->addInfluencedParticle(p);
        fBlackHole->addInfluencedParticle(p);
    }

    // update values from UI
    updateSimParams();

    // create the spatial hash
    spatialHash = new SpatialHash(2.1, system.getNumParticles()); // radius a little higher than 2*particle->radius
    spatialHash->create(system.getParticles());
}


void SceneFluid::updateSimParams()
{
    // get gravity from UI and update force
    double g = widget->getGravity();
    fGravity->setAcceleration(Vec3(0, -g, 0));

    // get other relevant UI values and update simulation params
    maxParticleLife = 20.0;
    emitRate = 100;
}


void SceneFluid::paint(const Camera& camera) {

    QOpenGLFunctions_3_3_Core* glFuncs = nullptr;
    glFuncs = QOpenGLContext::currentContext()->versionFunctions<QOpenGLFunctions_3_3_Core>();

    shader->bind();

    // camera matrices
    QMatrix4x4 camProj = camera.getPerspectiveMatrix();
    QMatrix4x4 camView = camera.getViewMatrix();
    shader->setUniformValue("ProjMatrix", camProj);
    shader->setUniformValue("ViewMatrix", camView);

    // lighting
    const int numLights = 1;
    const QVector3D lightPosWorld[numLights] = {QVector3D(100,500,100)};
    const QVector3D lightColor[numLights] = {QVector3D(1,1,1)};
    QVector3D lightPosCam[numLights];
    for (int i = 0; i < numLights; i++) {
        lightPosCam[i] = camView * lightPosWorld[i];
    }
    shader->setUniformValue("numLights", numLights);
    shader->setUniformValueArray("lightPos", lightPosCam, numLights);
    shader->setUniformValueArray("lightColor", lightColor, numLights);

    // Model matrix for transformations
    QMatrix4x4 modelMat;

    // draw the different spheres
    vaoSphereS->bind();

    // draw the water particles
    for (Particle* particle : waterParticles)
    {
        Vec3   p = particle->pos;
        Vec3   c = particle->color;
        double r = particle->radius;

        modelMat = QMatrix4x4();
        modelMat.translate(p[0], p[1], p[2]);
        modelMat.scale(r);
        shader->setUniformValue("ModelMatrix", modelMat);

        shader->setUniformValue("matdiff", GLfloat(c[0]), GLfloat(c[1]), GLfloat(c[2]));
        shader->setUniformValue("matshin", 100.0f);
        shader->setUniformValue("alpha", 1.0f);

        glFuncs->glDrawElements(GL_TRIANGLES, 3*numFacesSphereS, GL_UNSIGNED_INT, 0);
    }

    // Draw Floor
    vaoFloor->bind();
    modelMat.setToIdentity();
    modelMat.scale(40, 1, 40);
    shader->setUniformValue("ModelMatrix", modelMat);
    shader->setUniformValue("matdiff", 1.0f, 1.0f, 1.0f);
    shader->setUniformValue("alpha", 0.1f);
    glFuncs->glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

    // Draw Ceiling
    modelMat.setToIdentity();
    modelMat.translate(0, boundaryArea.y(), 0); // Translate to the top of the box
    modelMat.scale(40, 1, 40);
    shader->setUniformValue("ModelMatrix", modelMat);
    shader->setUniformValue("matdiff", 1.0f, 1.0f, 1.0f);
    shader->setUniformValue("alpha", 0.1f);
    glFuncs->glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

    // Draw Left Wall
    modelMat.setToIdentity();
    modelMat.translate(-boundaryArea.x(), boundaryArea.y() / 2, 0);
    modelMat.rotate(90, 0, 0, 1); // Rotate to stand vertically
    modelMat.scale(40, 1, 40);
    shader->setUniformValue("ModelMatrix", modelMat);
    shader->setUniformValue("matdiff", 1.0f, 1.0f, 1.0f);
    shader->setUniformValue("alpha", 0.1f);
    glFuncs->glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);


    // Draw Right Wall
    modelMat.setToIdentity();
    modelMat.translate(boundaryArea.x(), boundaryArea.y() / 2, 0);
    modelMat.rotate(90, 0, 0, 1); // Rotate to stand vertically
    modelMat.scale(40, 1, 40);
    shader->setUniformValue("ModelMatrix", modelMat);
    shader->setUniformValue("matdiff", 1.0f, 1.0f, 1.0f);
    shader->setUniformValue("alpha", 0.1f);
    glFuncs->glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

    // Draw Front Wall
    modelMat.setToIdentity();
    modelMat.translate(0, boundaryArea.y() / 2, boundaryArea.z());
    modelMat.rotate(90, 1, 0, 0); // Rotate to stand vertically
    modelMat.scale(40, 1, 40);
    shader->setUniformValue("ModelMatrix", modelMat);
    shader->setUniformValue("matdiff", 1.0f, 1.0f, 1.0f);
    shader->setUniformValue("alpha", 0.1f);
    glFuncs->glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

    // Draw Back Wall
    modelMat.setToIdentity();
    modelMat.translate(0, boundaryArea.y() / 2, -boundaryArea.z());
    modelMat.rotate(-90, 1, 0, 0); // Rotate to stand vertically
    modelMat.scale(40, 1, 40);
    shader->setUniformValue("ModelMatrix", modelMat);
    shader->setUniformValue("matdiff", 1.0f, 1.0f, 1.0f);
    shader->setUniformValue("alpha", 0.1f);
    glFuncs->glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

    // draw the boundary particles
//    for (Particle* particle : boundaryParticles)
//    {
//        Vec3   p = particle->pos;
//        Vec3   c = particle->color;
//        double r = particle->radius;

//        modelMat = QMatrix4x4();
//        modelMat.translate(p[0], p[1], p[2]);
//        modelMat.scale(r);
//        shader->setUniformValue("ModelMatrix", modelMat);

//        shader->setUniformValue("matdiff", 1.0f, 1.0f, 1.0f);
//        shader->setUniformValue("matspec", 0.1f, 0.1f, 0.1f);
//        shader->setUniformValue("matshin", 10.0f);
//        shader->setUniformValue("alpha", 0.05f);

//        glFuncs->glDrawElements(GL_TRIANGLES, 3*numFacesSphereS, GL_UNSIGNED_INT, 0);
//    }

    // draw the black hole
    if (blackHoleToggle)
    {
        vaoBlackHole->bind();
        modelMat = QMatrix4x4();
        modelMat.translate(blackHolePos[0], blackHolePos[1], blackHolePos[2]);
        modelMat.scale(blackHoleRadius);
        shader->setUniformValue("ModelMatrix", modelMat);

        shader->setUniformValue("matdiff", 0.0f, 0.0f, 0.0f);
        shader->setUniformValue("matspec", 1.0f, 1.0f, 1.0f);
        shader->setUniformValue("matshin", 100.f);

        glFuncs->glDrawElements(GL_TRIANGLES, 3*numFacesBlackHole, GL_UNSIGNED_INT, 0);
    }


    shader->release();
}

void SceneFluid::update(double dt) {

    // Update spatial hash with new positions
    spatialHash->create(system.getParticles());

    // 1. Find and store neighbors
    // Update neighbors using the spatialHash
    for (unsigned int i = 0; i < system.getNumParticles(); i++) {
        // Get the current particle, and clear its neighbors to begin
        Particle* particle = system.getParticle(i);
        particle->neighbors.clear();

        // Spatial hashing to find neighbors
        spatialHash->query(system.getParticles(), i, 2.1f);

        // Store the neighbors for the current particle
        for (unsigned int j = 0; j < spatialHash->querySize; j++) {
            int neighborIndex = spatialHash->queryIds[j];
            // Make sure the neighbor is not the current particle
            if (neighborIndex != i) {
                Particle* neighbor = system.getParticle(neighborIndex);
                particle->neighbors.push_back(neighbor);
            }
        }
    }

    // 2.-4.
    // SPH
    updateSPH();

    // Update the positions and velocities of the particles
    Vecd ppos = system.getPositions();
    integrator.step(system, dt);
    system.setPreviousPositions(ppos);

    // 5. Collisions
    // Water - Boundary
    for (unsigned int i = 0; i < system.getNumParticles(); i++) {
        Particle* particle = system.getParticle(i);

        for (Particle* neighbor : particle->neighbors) {
            // Perform collision test between particles
            if (testParticleCollision(particle, neighbor)) {
                resolveParticleCollision(particle, neighbor, kElastic, kFriction);
            }
        }

        // floor collision
        if (colliderFloor.testCollision(particle, dt)) {
            colliderFloor.resolveCollision(particle, kElastic, kFriction);
        }

        if (colliderCeiling.testCollision(particle, dt)) {
            colliderCeiling.resolveCollision(particle, kElastic, kFriction);
        }
        if (colliderLeftWall.testCollision(particle, dt)) {
            colliderLeftWall.resolveCollision(particle, kElastic, kFriction);
        }
        if (colliderRightWall.testCollision(particle, dt)) {
            colliderRightWall.resolveCollision(particle, kElastic, kFriction);
        }
        if (colliderFrontWall.testCollision(particle, dt)) {
            colliderFrontWall.resolveCollision(particle, kElastic, kFriction);
        }
        if (colliderBackWall.testCollision(particle, dt)) {
            colliderBackWall.resolveCollision(particle, kElastic, kFriction);
        }

    }

    // Update the system forces (must stay at the end of update)
    system.updateForces();
}




// -------------------------
// Helper functions for SPH

double getKernelPoly6(double r, double h){
    double h2_r2 = h*h - r*r;
    double h2_r2_cube = h2_r2*h2_r2*h2_r2;
    if(r >= 0 && r<=h){
        double h9 = h*h*h*h*h*h*h*h*h;
        return 315 / (64 * M_PI * h9) * h2_r2_cube;
    }
    return 0;
}


void SceneFluid::updateSPH() {
    // Calculate densities for each particle
    for (auto& p : system.getParticles()) {
        p->density = calculateDensity(p);
    }

    // Calculate forces and apply them
    for (auto& p : system.getParticles()) {
        Vec3 force = calculateForces(p);
        p->force += force;
    }
}

// 2. Density
double SceneFluid::calculateDensity(Particle* p) {
    double density = 0.0;
    for (Particle* neighbor : p->neighbors) {
        double distance = (neighbor->pos - p->pos).norm();
        if (distance < kernelRadius) {
            density += neighbor->mass * getKernelPoly6(distance, kernelRadius);
        }
    }
    p->density = density;
    return density;
}

// 3. Pressure
double SceneFluid::calculatePressure(double density) {
    return k_const * (density - restDensity);
}

// 4. Final SPH forces computations
Vec3 SceneFluid::calculateForces(Particle* p) {
    Vec3 totalForce(0.0, 0.0, 0.0);

    // Calculate pressure and viscosity forces
    for (Particle* neighbor : p->neighbors) {
        if (neighbor != p) {
            double distance = (neighbor->pos - p->pos).norm();
            if (distance < kernelRadius) {
                double densityPi = calculateDensity(p);
                double densityPj = calculateDensity(neighbor);
                double pressurePi = calculatePressure(densityPi) / (densityPi * densityPi);
                double pressurePj = calculatePressure(densityPj) / (densityPj * densityPj);

                Vec3 distanceVec = neighbor->pos - p->pos;

                // Calculate pressure force using the SPH pressure force formula
                Vec3 pressureForce = distanceVec.normalized() *
                                     (-neighbor->mass * (pressurePi + pressurePj) /
                                      (2.0 * densityPj) * getKernelPoly6(distance, kernelRadius));

                // Calculate viscosity force using the SPH viscosity force formula
                Vec3 viscosityForce = (neighbor->vel - p->vel) / densityPj *
                                      (2.0 * viscosityConstant * (kernelRadius - distance));

                // Add pressure and viscosity forces to the total force
                totalForce += pressureForce + viscosityForce;
            }
        }
    }

    return totalForce;
}

// end of Helper Functions
// -------------------------


void SceneFluid::mousePressed(const QMouseEvent* e, const Camera& cam)
{
    mouseX = e->pos().x();
    mouseY = e->pos().y();
}

void SceneFluid::mouseMoved(const QMouseEvent* e, const Camera& cam)
{
    int dx = e->pos().x() - mouseX;
    int dy = e->pos().y() - mouseY;
    mouseX = e->pos().x();
    mouseY = e->pos().y();

    Vec3 disp = cam.worldSpaceDisplacement(dx, -dy, cam.getEyeDistance());

    // example
    if (e->modifiers() & Qt::ControlModifier) {
        // move fountain
        // --removed--
    }
    else {
        Vec3 rayDir = cam.getRayDir(mouseX, mouseY);
        Vec3 rayOrigin = cam.getPos();

        // TODO: After I create my WCSPH, maybe handle some toys here again

//        bool hitAABB   = checkIntersectionWithAABB(rayOrigin, rayDir, colliderPropBox);
//        bool hitSphere = checkIntersectionWithSphere(rayOrigin, rayDir, colliderPropSphere);

//        if (hitAABB)
//            colliderPropBox.Move(disp);
//        if (hitSphere)
//            colliderPropSphere.Move(disp);
    }
}


bool SceneFluid::checkIntersectionWithAABB(const Vec3& rayOrigin, const Vec3& rayDir, const ColliderAABB& collider)
{
    Eigen::Array3d tMinArray = (collider.aabbPosition - collider.aabbScale - rayOrigin).array() / rayDir.array();
    Eigen::Array3d tMaxArray = (collider.aabbPosition + collider.aabbScale - rayOrigin).array() / rayDir.array();

    Eigen::Array3d t1Array = (tMinArray < tMaxArray).select(tMinArray, tMaxArray);  // min of the 2
    Eigen::Array3d t2Array = (tMinArray < tMaxArray).select(tMaxArray, tMinArray);  // max of the 2


    float tNear = std::max(std::max(t1Array.x(), t1Array.y()), t1Array.z());    // max of the 3
    float tFar =  std::min(std::min(t2Array.x(), t2Array.y()), t2Array.z());    // min of the 3

    if (tNear <= tFar)
        return true;

    return false;
}

bool SceneFluid::checkIntersectionWithSphere(const Vec3& rayOrigin, const Vec3& rayDir, const ColliderSphere& collider)
{
    Vec3 radius = Vec3(collider.radius, collider.radius, collider.radius);
    Eigen::Array3d tMinArray = (collider.center - radius - rayOrigin).array() / rayDir.array();
    Eigen::Array3d tMaxArray = (collider.center + radius - rayOrigin).array() / rayDir.array();

    Eigen::Array3d t1Array = (tMinArray < tMaxArray).select(tMinArray, tMaxArray);  // min of the 2
    Eigen::Array3d t2Array = (tMinArray < tMaxArray).select(tMaxArray, tMinArray);  // max of the 2


    float tNear = std::max(std::max(t1Array.x(), t1Array.y()), t1Array.z());    // max of the 3
    float tFar =  std::min(std::min(t2Array.x(), t2Array.y()), t2Array.z());    // min of the 3

    if (tNear <= tFar)
        return true;

    return false;
}

// Water particle creation

void SceneFluid::createWaterParticles() {
    const float waterParticleSpacing = 2.0f; // Spacing between water particles
    const float particleMass = 0.01f;
    const Vec3 blue(0.f, 0.f, 1.f);

    // Start and End of the Water Dam
    // Essentially handling the number of water particles
    Vec3 waterDamStart(-boundaryArea.x() + 0.5f, 2.f, -boundaryArea.z());
    Vec3 waterDamEnd(-boundaryArea.x()/2, 25.f, -boundaryArea.z()/2);

    // Create water particles within the water dam dimensions
    for (float x = waterDamStart.x(); x <= waterDamEnd.x(); x += waterParticleSpacing) {
        for (float y = waterDamStart.y(); y <= waterDamEnd.y(); y += waterParticleSpacing) {
            for (float z = waterDamStart.z(); z <= waterDamEnd.z(); z += waterParticleSpacing) {
                Vec3 position(x, y, z);
                Particle* waterParticle = new Particle(position);
                waterParticle->color = blue;
                waterParticle->mass = particleMass;
                waterParticle->boundary = false;
                waterParticle->density = restDensity;
                waterParticle->pressure = 0.0f; // Initial pressure is zero
                waterParticles.push_back(waterParticle);
                system.addParticle(waterParticle);

                fGravity->addInfluencedParticle(waterParticle);
                fBlackHole->addInfluencedParticle(waterParticle);
            }
        }
    }

    // Some additional balls above the previous to create some ripple effect
    int numLayers = 5;
    int particlesPerLayer = 5;
    // Create additional particles above the existing layers
    for (int layer = 0; layer < numLayers; layer++) {
        for (int i = 0; i < particlesPerLayer; i++) {
            // Randomly position the particles above the existing layers
            float x = waterDamStart.x() + static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / (waterDamEnd.x() - waterDamStart.x())));
            float y = waterDamEnd.y() + (layer + 1) * waterParticleSpacing; // Position above existing layers
            float z = waterDamStart.z() + static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / (waterDamEnd.z() - waterDamStart.z())));

            Vec3 position(x, y, z);
            Particle* waterParticle = new Particle(position);
            waterParticle->color = blue;
            waterParticle->mass = particleMass;
            waterParticle->boundary = false;
            waterParticle->density = restDensity;
            waterParticle->pressure = 0.0f; // Initial pressure is zero
            waterParticles.push_back(waterParticle);
            system.addParticle(waterParticle);

            fGravity->addInfluencedParticle(waterParticle);
            fBlackHole->addInfluencedParticle(waterParticle);
        }
    }
}


/*
* NOT USED AFTER ALL, I USED 6 PLANES
*/
// Boundary particle creation
void SceneFluid::createBoundaryParticles()
{
    const float boundaryParticleSpacing = 1.0f;
    const float particleMass = 1.0f;

    // Calculate start and end positions for particles in each dimension
    Vec3 startPos = colliderPropBox.aabbPosition - colliderPropBox.aabbScale;
    Vec3 endPos = colliderPropBox.aabbPosition + colliderPropBox.aabbScale;

    // Iterate over each face of the AABB and create particles
    for (float x = startPos.x(); x <= endPos.x(); x += boundaryParticleSpacing) {
        for (float y = startPos.y(); y <= endPos.y(); y += boundaryParticleSpacing) {
            for (float z = startPos.z(); z <= endPos.z(); z += boundaryParticleSpacing) {
                // Skip internal points, only place particles on the faces
                if (x > startPos.x() && x < endPos.x() && y > startPos.y() && y < endPos.y() && z > startPos.z() && z < endPos.z()) {
                    continue;
                }

                Vec3 position(x, y, z);
                Particle* boundaryParticle = new Particle(position);
                boundaryParticle->color = Vec3(0.0f, 1.0f, 1.0f);
                boundaryParticle->mass = particleMass;
                boundaryParticle->boundary = true;
                boundaryParticle->density = restDensity;
                boundaryParticle->pressure = 0.0f; // Initial pressure is zero
                boundaryParticles.push_back(boundaryParticle);
                system.addParticle(boundaryParticle);
            }
        }
    }
}
