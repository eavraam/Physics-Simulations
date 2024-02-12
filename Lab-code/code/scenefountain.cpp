#include "scenefountain.h"
#include "glutils.h"
#include "model.h"
#include <QOpenGLFunctions_3_3_Core>


SceneFountain::SceneFountain() {
    widget = new WidgetFountain();
    connect(widget, SIGNAL(updatedParameters()), this, SLOT(updateSimParams()));
}


SceneFountain::~SceneFountain() {
    if (widget)         delete widget;
    if (shader)         delete shader;
    if (fGravity)       delete fGravity;
    if (vaoFloor)       delete vaoFloor;
    if (vaoSphereS)     delete vaoSphereS;
    if (vaoPropSphere)  delete vaoPropSphere;
    if (vaoPropBox)     delete vaoPropBox;
    if (vaoBlackHole)   delete vaoBlackHole;
}


void SceneFountain::initialize(double dt, int dr, double kEl, double kFr, int bht)
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

    // scene
    fountainPos = Vec3(0, 5, 0);
    colliderFloor.setPlane(Vec3(0, 1, 0), 0);
    colliderPropSphere.setSphere(Vec3(-30.0f, 0.0f, -10.0f), 20.0f);
    colliderPropBox.setAABB(Vec3(15.0f, 50.0f, 0.0f), Vec3(20.0f, 10.0f, 10.0f));

    // store all the PROP colliders in a vector
    colliders.push_back(&colliderPropBox);
    colliders.push_back(&colliderPropSphere);

    // create the spatial hash
    //spatialHash = new SpatialHash(2, 2000);
}


void SceneFountain::reset(double dt, int dr, double kEl, double kFr, int bht)
{

    // Reset drag, elasticity, friction and black hole
    timeStep = dt;
    drag = dr;
    kElastic = kEl;
    kFriction = kFr;
    blackHoleToggle = bht;

    // reset the spatial hash
    //spatialHash->create(system.getParticles());

    std::cout << "Reset kElastic to: " << kElastic << std::endl;
    std::cout << "Reset kFriction to: " << kFriction << std::endl;
    std::cout << "Reset drag to: " << drag << std::endl;
    std::cout << "Reset black hole toggle to: " << blackHoleToggle << std::endl;

    // Forces
    // reset system forces
    system.clearForces();


    // gravity
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

    // update values from UI
    updateSimParams();

    // reset random seed
    Random::seed(1337);

    // erase all particles
    fGravity->clearInfluencedParticles();
    // NULL pointer reference, need to handle it.
    if (fDragLinear != nullptr)
        fDragLinear->clearInfluencedParticles();
    if (fDragQuadratic != nullptr)
        fDragQuadratic->clearInfluencedParticles();
    if (fBlackHole != nullptr)
        fBlackHole->clearInfluencedParticles();
    system.deleteParticles();
    deadParticles.clear();
}


void SceneFountain::updateSimParams()
{
    // get gravity from UI and update force
    double g = widget->getGravity();
    fGravity->setAcceleration(Vec3(0, -g, 0));

    // get other relevant UI values and update simulation params
    maxParticleLife = 20.0;
    emitRate = 10;
}


void SceneFountain::paint(const Camera& camera) {

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

    // draw floor
    vaoFloor->bind();
    QMatrix4x4 modelMat;
    modelMat.scale(100, 1, 100);
    shader->setUniformValue("ModelMatrix", modelMat);
    shader->setUniformValue("matdiff", 0.8f, 0.8f, 0.8f);
    shader->setUniformValue("matspec", 0.0f, 0.0f, 0.0f);
    shader->setUniformValue("matshin", 0.0f);
    glFuncs->glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

    // draw the different spheres
    vaoSphereS->bind();
    for (const Particle* particle : system.getParticles()) {
        Vec3   p = particle->pos;
        Vec3   c = particle->color;
        double r = particle->radius;

        modelMat = QMatrix4x4();
        modelMat.translate(p[0], p[1], p[2]);
        modelMat.scale(r);
        shader->setUniformValue("ModelMatrix", modelMat);

        shader->setUniformValue("matdiff", GLfloat(c[0]), GLfloat(c[1]), GLfloat(c[2]));
        shader->setUniformValue("matspec", 1.0f, 1.0f, 1.0f);
        shader->setUniformValue("matshin", 100.f);

        glFuncs->glDrawElements(GL_TRIANGLES, 3*numFacesSphereS, GL_UNSIGNED_INT, 0);
    }

    // draw the different props
    // Sphere prop
    vaoPropSphere->bind(); // already bound, but for order-safety reasons
    modelMat = QMatrix4x4();
    modelMat.translate(colliderPropSphere.center.x(),
                       colliderPropSphere.center.y(),
                       colliderPropSphere.center.z()
                       );
    modelMat.scale(colliderPropSphere.radius);
    shader->setUniformValue("ModelMatrix", modelMat);

    shader->setUniformValue("matdiff", 1.0f, 0.0f, 0.0f);
    shader->setUniformValue("matspec", 1.0f, 1.0f, 1.0f);
    shader->setUniformValue("matshin", 100.f);

    glFuncs->glDrawElements(GL_TRIANGLES, 3*numFacesPropSphere, GL_UNSIGNED_INT, 0);

    // Box prop
    vaoPropBox->bind();
    modelMat = QMatrix4x4();
    modelMat.translate(colliderPropBox.aabbPosition.x(),
                       colliderPropBox.aabbPosition.y(),
                       colliderPropBox.aabbPosition.z()
                       );
    modelMat.scale(
                colliderPropBox.aabbScale.x(),
                colliderPropBox.aabbScale.y(),
                colliderPropBox.aabbScale.z()
                );
    shader->setUniformValue("ModelMatrix", modelMat);

    shader->setUniformValue("matdiff", 0.65f, 0.65f, 0.65f);
    shader->setUniformValue("matspec", 1.0f, 1.0f, 1.0f);
    shader->setUniformValue("matshin", 100.f);

    glFuncs->glDrawElements(GL_TRIANGLES, 3*numFacesPropBox, GL_UNSIGNED_INT, 0);


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


void SceneFountain::update(double dt) {

    // emit new particles, reuse dead ones if possible
    int emitParticles = std::max(1, int(std::round(emitRate * dt)));
    for (int i = 0; i < emitParticles; i++) {
        Particle* p;
        if (!deadParticles.empty()) {
            // reuse one dead particle
            p = deadParticles.front();
            deadParticles.pop_front();
        }
        else {
            // create new particle
            p = new Particle();
            system.addParticle(p);

            // don't forget to add particle to forces that affect it
            fGravity->addInfluencedParticle(p);

            if (blackHoleToggle)
            {
                fBlackHole->addInfluencedParticle(p);
            }
            if (drag == 1)
            {
                fDragLinear->addInfluencedParticle(p);
            }
            else if (drag == 2)
            {
                fDragQuadratic->addInfluencedParticle(p);
            }
        }

        p->color = Vec3(153/255.0, 217/255.0, 234/255.0);
        p->radius = 1.0;
        p->life = maxParticleLife;

        double x = 0;
        double y = 0;
        double z = 0;
        p->pos = Vec3(x, y, z) + fountainPos;
        p->vel = Vec3(
                    Random::get(0.0, 10.0) - 5,
                    30,
                    Random::get(0.0, 10.0) - 5
                );
    }

    // create the spatial hash
    //spatialHash->create(system.getParticles());

    // integration step
    Vecd ppos = system.getPositions();
    integrator.step(system, dt);
    system.setPreviousPositions(ppos);

    // collisions
    for (Particle* p : system.getParticles()) {
        // floor collision
        if (colliderFloor.testCollision(p, dt)) {
            colliderFloor.resolveCollision(p, kElastic, kFriction);
        }
        // sphere prop collision
        if (colliderPropSphere.testCollision(p, dt)) {
            colliderPropSphere.resolveCollision(p, kElastic, kFriction);
        }
        // box prop collision
        if (colliderPropBox.testCollision(p, dt)) {
           colliderPropBox.resolveCollision(p, kElastic, kFriction);
        }

//        // spatial hash collision
//        float particleMinDist = 2.0f * p->radius;

//        // query the spatial hash for collisions
//        spatialHash->query(system.getParticles(), p->id, particleMinDist);

//        // interball collision
//        // particle-particle collision using spatial hashing
//        for(unsigned int nr=0; nr<spatialHash->querySize;nr++){
//            Particle *p2 = system.getParticles()[spatialHash->queryIds[nr]];
//            Vecd normal = p->pos - p2->pos;
//            double d  = (normal).norm();
//            double d2 = d*d;

//            // are the balls overlapping?
//            if(d2 > 0.f && d2 < particleMinDist*particleMinDist) {
//                normal = normal/d;

//                double corr = (particleMinDist - d) * 0.5;

//                // separate the balls
//                p->pos  += normal*corr;
//                p2->pos -= normal*corr;

//                // reflect velocities along normal
//                double vi = p->vel.dot(normal);
//                double vj = p2->vel.dot(normal);

//                p->vel  += normal*(vj-vi);
//                p2->vel += normal*(vi-vj);

//                p->color  = Vec3(1, 0, 0);
//                p2->color = Vec3(1, 0, 0);

//            }
//        }
//        p->color = Vec3(153/255.0, 217/255.0, 234/255.0);

    }

    // check dead particles
    for (Particle* p : system.getParticles()) {
        if (p->life > 0) {
            p->life -= dt;
            if (p->life < 0) {
                deadParticles.push_back(p);
            }
        }
    }
}

void SceneFountain::mousePressed(const QMouseEvent* e, const Camera& cam)
{
    mouseX = e->pos().x();
    mouseY = e->pos().y();
}

void SceneFountain::mouseMoved(const QMouseEvent* e, const Camera& cam)
{
    int dx = e->pos().x() - mouseX;
    int dy = e->pos().y() - mouseY;
    mouseX = e->pos().x();
    mouseY = e->pos().y();

    Vec3 disp = cam.worldSpaceDisplacement(dx, -dy, cam.getEyeDistance());

    // example
    if (e->modifiers() & Qt::ControlModifier) {
        // move fountain
        fountainPos += disp;
    }
    else {
        Vec3 rayDir = cam.getRayDir(mouseX, mouseY);
        Vec3 rayOrigin = cam.getPos();

        bool hitAABB   = checkIntersectionWithAABB(rayOrigin, rayDir, colliderPropBox);
        bool hitSphere = checkIntersectionWithSphere(rayOrigin, rayDir, colliderPropSphere);

        if (hitAABB)
            colliderPropBox.Move(disp);
        if (hitSphere)
            colliderPropSphere.Move(disp);
    }
}


bool SceneFountain::checkIntersectionWithAABB(const Vec3& rayOrigin, const Vec3& rayDir, const ColliderAABB& collider)
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

bool SceneFountain::checkIntersectionWithSphere(const Vec3& rayOrigin, const Vec3& rayDir, const ColliderSphere& collider)
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

