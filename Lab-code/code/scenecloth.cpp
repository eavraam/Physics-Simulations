#include "scenecloth.h"
#include "glutils.h"
#include "model.h"
#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLBuffer>


SceneCloth::SceneCloth() {
    widget = new WidgetCloth();
    connect(widget, SIGNAL(updatedParameters()), this, SLOT(updateSimParams()));
    connect(widget, SIGNAL(freeAnchors()), this, SLOT(freeAnchors()));
}

SceneCloth::~SceneCloth() {
    if (widget)       delete widget;
    if (shaderPhong)  delete shaderPhong;
    if (vaoSphereS)   delete vaoSphereS;
    if (vaoSphereL)   delete vaoSphereL;
    if (vaoFloor)     delete vaoFloor;
    if (vaoCube)      delete vaoCube;
    if (vaoMesh)      delete vaoMesh;
    if (vaoBlackHole) delete vaoBlackHole;
    if (vboMesh)      delete vboMesh;
    if (iboMesh)      delete iboMesh;

    system.deleteParticles();
    if (fGravity)     delete fGravity;
    if (fBlackHole)   delete fBlackHole;
}

void SceneCloth::initialize(double dt, int dr, double kEl, double kFr, int bht) {

    // Initialize drag, elasticity, friction and black hole
    timeStep = dt;
    drag = dr;
    kElastic = kEl;
    kFriction = kFr;
    blackHoleToggle = bht;

    // load shaders
    shaderPhong = glutils::loadShaderProgram(":/shaders/phong.vert", ":/shaders/phong.frag");
    shaderCloth = glutils::loadShaderProgram(":/shaders/cloth.vert", ":/shaders/cloth.geom", ":/shaders/cloth.frag");

    // create floor VAO
    Model quad = Model::createQuad();
    vaoFloor = glutils::createVAO(shaderPhong, &quad);
    glutils::checkGLError();

    // create sphere VAOs
    Model sphere = Model::createIcosphere(3);
    vaoSphereL = glutils::createVAO(shaderPhong, &sphere);
    numFacesSphereL = sphere.numFaces();
    glutils::checkGLError();

    sphere = Model::createIcosphere(1);
    vaoSphereS = glutils::createVAO(shaderPhong, &sphere);
    numFacesSphereS = sphere.numFaces();
    glutils::checkGLError();

    // create props' VAOs
    Model propSphere = Model::createIcosphere(1);
    vaoPropSphere = glutils::createVAO(shaderPhong, &propSphere);
    numFacesPropSphere = propSphere.numFaces();
    glutils::checkGLError();

    Model propBox = Model::createCube();
    vaoPropBox = glutils::createVAO(shaderPhong, &propBox);
    numFacesPropBox = propBox.numFaces();
    glutils::checkGLError();

    // create cube VAO
    Model cube = Model::createCube();
    vaoCube = glutils::createVAO(shaderPhong, &cube);
    glutils::checkGLError();

    // create black hole VAO
    Model blackHole = Model::createIcosphere(1);
    vaoBlackHole = glutils::createVAO(shaderPhong, &blackHole);
    numFacesBlackHole = blackHole.numFaces();
    glutils::checkGLError();

    // create gravity force
    fGravity = new ForceConstAcceleration();
    system.addForce(fGravity);

    // black hole
    blackHolePos = Vec3(-80.0f, 100.0f, -10.0f);
    blackHoleRadius = 5.0f;
    blackHoleMass = 10000.0f;

    fBlackHole = new ForceBlackHole(blackHolePos, blackHoleMass);
    if (blackHoleToggle)
    {
        system.addForce(fBlackHole);
    }

    // create cloth mesh VAO
    vaoMesh = new QOpenGLVertexArrayObject();
    vaoMesh->create();
    vaoMesh->bind();
    vboMesh = new QOpenGLBuffer(QOpenGLBuffer::Type::VertexBuffer);
    vboMesh->create();
    vboMesh->bind();
    vboMesh->setUsagePattern(QOpenGLBuffer::UsagePattern::DynamicDraw);
    vboMesh->allocate(1000*1000*3*3*sizeof(float)); // sync with widget max particles
    shaderCloth->setAttributeBuffer("vertex", GL_FLOAT, 0, 3, 0);
    shaderCloth->enableAttributeArray("vertex");
    iboMesh = new QOpenGLBuffer(QOpenGLBuffer::Type::IndexBuffer);
    iboMesh->create();
    iboMesh->bind();
    iboMesh->setUsagePattern(QOpenGLBuffer::UsagePattern::StaticDraw);
    iboMesh->allocate(1000*1000*2*3*sizeof(unsigned int));
    vaoMesh->release();
}

void SceneCloth::reset(double dt, int dr, double kEl, double kFr, int bht)
{

    // Reset drag, elasticity, friction and black hole
    timeStep = dt;
    drag = dr;
    kElastic = kEl;
    kFriction = kFr;
    blackHoleToggle = bht;

    // we only update numParticles on resets
    updateSimParams();

    // reset forces
    fGravity->clearInfluencedParticles();
    fBlackHole->clearInfluencedParticles();
    for(int i=0; i<sceneForceSprings.length(); i++)
    {
        sceneForceSprings[i]->clearInfluencedParticles();
    }
    system.clearForces();
    sceneForceSprings.clear();

    // reset particles
    system.deleteParticles();
    sceneParticles.clear();

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

    // cloth params
    Vec2 dims = widget->getDimensions();
    Vec2i dimParticles = widget->getNumParticles();
    numParticlesX = dimParticles.x();
    numParticlesY = dimParticles.y();
    particleRadius = widget->getParticleRadius();
    // compute cloth
    computeCloth(numParticlesX, numParticlesY, sceneParticles, sceneForceSprings, dims[0], dims[1], particleRadius);

    // create particles
    numParticles = numParticlesX * numParticlesY;
    fixedParticle = std::vector<bool>(numParticles, false);

    for (int i=0; i < sceneParticles.length(); i++)
    {
        if (sceneParticles[i]->isFixed())
        {
            fixedParticle[i] = true;
        }
        system.addParticle(sceneParticles[i]);
        fGravity->addInfluencedParticle(sceneParticles[i]);
        fBlackHole->addInfluencedParticle(sceneParticles[i]);
    }

    // Add spring forces
    for (ForceSpring* fs : sceneForceSprings) {
        system.addForce(fs);
    }

    // Code for PROVOT layout
    updateSprings();

    // update index buffer
    iboMesh->bind();
    numMeshIndices = (numParticlesX - 1)*(numParticlesY - 1)*2*3;
    int* indices = new int[numMeshIndices];
    int idx = 0;
    for (int i = 0; i < numParticlesX-1; i++) {
        for (int j = 0; j < numParticlesY-1; j++) {
            indices[idx  ] = i*numParticlesY + j;
            indices[idx+1] = (i+1)*numParticlesY + j;
            indices[idx+2] = i*numParticlesY + j + 1;
            indices[idx+3] = i*numParticlesY + j + 1;
            indices[idx+4] = (i+1)*numParticlesY + j;
            indices[idx+5] = (i+1)*numParticlesY + j + 1;
            idx += 6;
        }
    }
    void* bufptr = iboMesh->mapRange(0, numMeshIndices*sizeof(int),
                                     QOpenGLBuffer::RangeInvalidateBuffer | QOpenGLBuffer::RangeWrite);
    memcpy(bufptr, (void*)(indices), numMeshIndices*sizeof(int));
    iboMesh->unmap();
    iboMesh->release();
    delete[] indices;
    glutils::checkGLError();

    // Add collider prop objects
    colliderFloor.setPlane(Vec3(0, 1, 0), 0);
    colliderPropSphere.setSphere(Vec3(0.0f, 27.0f, -10.0f), 20.0f);
//    colliderPropBox.setAABB(Vec3(30.0f, 0.0f, 0.0f), Vec3(20.0f, 10.0f, 10.0f));

    // Create spatial hash
    spatial_hash_toggle = widget->useSpatialHash();
    spatialHash = new SpatialHash(2.5f, numParticles); // spacing a little more than 2*radius
    spatialHash->create(system.getParticles());
}


void SceneCloth::updateSprings()
{
    double ks = widget->getStiffness();
    double kd = widget->getDamping();

    // here I update all ks and kd parameters.
    // idea: if you want to enable/disable a spring type, you can set ks to 0 for these
    for (ForceSpring* f : sceneForceSprings)
    {
        f->ks = ks;
        f->kd = kd;
    }
}

void SceneCloth::updateSimParams()
{
    double g = widget->getGravity();
    fGravity->setAcceleration(Vec3(0, -g, 0));

    updateSprings();

    for (Particle* p : system.getParticles()) {
        p->radius = widget->getParticleRadius();
    }

    showParticles = widget->showParticles();
    spatial_hash_toggle = widget->useSpatialHash();
}

void SceneCloth::freeAnchors()
{
    fixedParticle = std::vector<bool>(numParticles, false);
}

void SceneCloth::paint(const Camera& camera)
{
    QOpenGLFunctions_3_3_Core* glFuncs = nullptr;
    glFuncs = QOpenGLContext::currentContext()->versionFunctions<QOpenGLFunctions_3_3_Core>();

    shaderPhong->bind();
    shaderPhong->setUniformValue("normalSign", 1.0f);

    // camera matrices
    QMatrix4x4 camProj = camera.getPerspectiveMatrix();
    QMatrix4x4 camView = camera.getViewMatrix();
    shaderPhong->setUniformValue("ProjMatrix", camProj);
    shaderPhong->setUniformValue("ViewMatrix", camView);

    // lighting
    const int numLights = 1;
    const QVector3D lightPosWorld[numLights] = {QVector3D(80,80,80)};
    const QVector3D lightColor[numLights] = {QVector3D(1,1,1)};
    QVector3D lightPosCam[numLights];
    for (int i = 0; i < numLights; i++) {
        lightPosCam[i] = camView * lightPosWorld[i];
    }
    shaderPhong->setUniformValue("numLights", numLights);
    shaderPhong->setUniformValueArray("lightPos", lightPosCam, numLights);
    shaderPhong->setUniformValueArray("lightColor", lightColor, numLights);

    // Model matrix
    QMatrix4x4 modelMat;

    // draw floor
    vaoFloor->bind();
    modelMat.scale(100, 1, 100);
    shaderPhong->setUniformValue("ModelMatrix", modelMat);
    shaderPhong->setUniformValue("matdiff", 0.8f, 0.8f, 0.8f);
    shaderPhong->setUniformValue("matspec", 0.0f, 0.0f, 0.0f);
    shaderPhong->setUniformValue("matshin", 0.0f);
    glFuncs->glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

    // draw the particle spheres
    if (showParticles) {
        vaoSphereS->bind();
        shaderPhong->setUniformValue("matspec", 1.0f, 1.0f, 1.0f);
        shaderPhong->setUniformValue("matshin", 100.f);
        for (int i = 0; i < numParticles; i++) {
            const Particle* particle = system.getParticle(i);
            Vec3   p = particle->pos;
            Vec3   c = particle->color;
            if (fixedParticle[i])      c = Vec3(63/255.0, 72/255.0, 204/255.0);
            if (i == selectedParticle) c = Vec3(1.0,0.9,0);

            modelMat = QMatrix4x4();
            modelMat.translate(p[0], p[1], p[2]);
            modelMat.scale(particle->radius);
            shaderPhong->setUniformValue("ModelMatrix", modelMat);
            shaderPhong->setUniformValue("matdiff", GLfloat(c[0]), GLfloat(c[1]), GLfloat(c[2]));
            glFuncs->glDrawElements(GL_TRIANGLES, 3*numFacesSphereS, GL_UNSIGNED_INT, 0);
        }
    }

    // TODO: draw colliders and walls

    // draw the different props
    // Sphere prop
    vaoPropSphere->bind(); // already bound, but for order-safety reasons
    modelMat = QMatrix4x4();
    modelMat.translate(colliderPropSphere.center.x(),
                       colliderPropSphere.center.y(),
                       colliderPropSphere.center.z()
                       );
    modelMat.scale(colliderPropSphere.radius);
    shaderPhong->setUniformValue("ModelMatrix", modelMat);

    shaderPhong->setUniformValue("matdiff", 0.65f, 0.65f, 0.65f);
    shaderPhong->setUniformValue("matspec", 1.0f, 1.0f, 1.0f);
    shaderPhong->setUniformValue("matshin", 100.f);

    glFuncs->glDrawElements(GL_TRIANGLES, 3*numFacesPropSphere, GL_UNSIGNED_INT, 0);

//    // Box prop
//    vaoPropBox->bind();
//    modelMat = QMatrix4x4();
//    modelMat.translate(colliderPropBox.aabbPosition.x(),
//                       colliderPropBox.aabbPosition.y(),
//                       colliderPropBox.aabbPosition.z()
//                       );
//    modelMat.scale(
//                colliderPropBox.aabbScale.x(),
//                colliderPropBox.aabbScale.y(),
//                colliderPropBox.aabbScale.z()
//                );
//    shaderPhong->setUniformValue("ModelMatrix", modelMat);

//    shaderPhong->setUniformValue("matdiff", 0.65f, 0.65f, 0.65f);
//    shaderPhong->setUniformValue("matspec", 1.0f, 1.0f, 1.0f);
//    shaderPhong->setUniformValue("matshin", 100.f);

//    glFuncs->glDrawElements(GL_TRIANGLES, 3*numFacesPropBox, GL_UNSIGNED_INT, 0);

    // draw the black hole
    if (blackHoleToggle)
    {
        vaoBlackHole->bind();
        modelMat = QMatrix4x4();
        modelMat.translate(blackHolePos[0], blackHolePos[1], blackHolePos[2]);
        modelMat.scale(blackHoleRadius);
        shaderPhong->setUniformValue("ModelMatrix", modelMat);

        shaderPhong->setUniformValue("matdiff", 0.0f, 0.0f, 0.0f);
        shaderPhong->setUniformValue("matspec", 1.0f, 1.0f, 1.0f);
        shaderPhong->setUniformValue("matshin", 100.f);

        glFuncs->glDrawElements(GL_TRIANGLES, 3*numFacesBlackHole, GL_UNSIGNED_INT, 0);
    }


    // Release shader
    shaderPhong->release();


    // update cloth mesh VBO coords
    vboMesh->bind();
    float* pos = new float[3*numParticles];
    for (int i = 0; i < numParticles; i++) {
        pos[3*i  ] = system.getParticle(i)->pos.x();
        pos[3*i+1] = system.getParticle(i)->pos.y();
        pos[3*i+2] = system.getParticle(i)->pos.z();
    }
    void* bufptr = vboMesh->mapRange(0, 3*numParticles*sizeof(float),
                       QOpenGLBuffer::RangeInvalidateBuffer | QOpenGLBuffer::RangeWrite);
    memcpy(bufptr, (void*)(pos), 3*numParticles*sizeof(float));
    vboMesh->unmap();
    vboMesh->release();
    delete[] pos;

    // draw mesh
    shaderCloth->bind();
    shaderCloth->setUniformValue("ProjMatrix", camProj);
    shaderCloth->setUniformValue("ViewMatrix", camView);
    shaderCloth->setUniformValue("NormalMatrix", camView.normalMatrix());
    shaderCloth->setUniformValue("matdiffFront", 0.7f, 0.0f, 0.0f);
    shaderCloth->setUniformValue("matspecFront", 1.0f, 1.0f, 1.0f);
    shaderCloth->setUniformValue("matshinFront", 100.0f);
    shaderCloth->setUniformValue("matdiffBack", 0.7f, 0.3f, 0.0f);
    shaderCloth->setUniformValue("matspecBack", 0.0f, 0.0f, 0.0f);
    shaderCloth->setUniformValue("matshinBack", 0.0f);
    shaderCloth->setUniformValue("numLights", numLights);
    shaderCloth->setUniformValueArray("lightPos", lightPosCam, numLights);
    shaderCloth->setUniformValueArray("lightColor", lightColor, numLights);
    vaoMesh->bind();
    glFuncs->glDrawElements(GL_TRIANGLES, numMeshIndices, GL_UNSIGNED_INT, 0);
    vaoMesh->release();
    shaderCloth->release();

    glutils::checkGLError();
}


void SceneCloth::update(double dt)
{
    // Create spatial hash
    spatialHash->create(system.getParticles());

    // fixed particles: no velocity, no force acting
    for (int i = 0; i < numParticles; i++) {
        if (fixedParticle[i]) {
            Particle* p = system.getParticle(i);
            p->vel = Vec3(0,0,0);
            p->force = Vec3(0,0,0);
        }
    }

    // integration step
    Vecd ppos = system.getPositions();
    integrator.step(system, dt);
    system.setPreviousPositions(ppos);

    // user interaction
    // TODO: test and resolve for collisions during user movement
    if (selectedParticle >= 0) {
        Particle* p = system.getParticle(selectedParticle);
        p->pos  = cursorWorldPos;
        p->vel = Vec3(0,0,0);
    }

    // TODO: relaxation
    // Relaxation
    int counter = 0;
    while (counter < relaxationSteps) {
        for (ForceSpring* fs : sceneForceSprings) {

            std::vector<Particle*> fs_particles = fs->getInfluencedParticles();
            Vec3 direction = fs_particles[1]->pos - fs_particles[0]->pos;
            float currentLength = direction.norm();
            float dirPrev = fs->L_distance;

            // Calculate the adjustment ratio (slight stretch or compress)
            float adjustRatio;
            if (currentLength > 1.1 * dirPrev) {
                adjustRatio = 1.05; // Slight stretching
            } else if (currentLength < 0.9 * dirPrev) {
                adjustRatio = 0.95; // Slight compression
            } else {
                // No adjustment needed
                continue;
            }

            Vec3 adjustedDir = direction.normalized() * (dirPrev * adjustRatio);
            Vec3 center = (fs_particles[1]->pos + fs_particles[0]->pos) / 2.0f;

            // Adjust the positions of the particles
            if (!fs_particles[0]->isFixed()) {
                fs_particles[0]->pos = center - adjustedDir / 2.0f;
            }
            if (!fs_particles[1]->isFixed()) {
                fs_particles[1]->pos = center + adjustedDir / 2.0f;
            }
        }
        counter++;
    }

    // TODO: test and resolve collisions
    // collisions
    for (int i=0; i<numParticles; i++) {

        Particle* p = system.getParticle(i);

        // Spatial hashing to prevent self-intersections
        if (spatial_hash_toggle)
        {
            // Hash testing
            spatialHash->query(system.getParticles(), i, collisionRadius);

            for (unsigned int j = 0; j < spatialHash->querySize; j++) {
                int neighborIdx = spatialHash->queryIds[j];
                if (i == neighborIdx) continue; // Skip self

                Particle* neighbor = system.getParticle(neighborIdx);

                // Check and resolve collision
                if (testParticleCollision(p, neighbor)) {
                    resolveParticleCollision(p, neighbor, kElastic, kFriction);
                }
            }
        }

        // floor collision
        if (colliderFloor.testCollision(p, dt)) {
            colliderFloor.resolveCollision(p, kElastic, kFriction);
        }
        // sphere prop collision
        if (colliderPropSphere.testCollision(p, dt)) {
            colliderPropSphere.resolveCollision(p, kElastic, kFriction);
        }
//        // box prop collision
//        if (colliderPropBox.testCollision(p, dt)) {
//           colliderPropBox.resolveCollision(p, kElastic, kFriction);
//        }
    }

    // Update velocities, since the damping term of Spring depends on the velocity
    for (int i = 0; i < numParticles; i++) {
        Particle* p = system.getParticle(i);
        p->vel = (p->pos - p->prevPos) / dt;
    }

    // needed after we have done collisions and relaxation, since spring forces depend on p and v
    system.updateForces();
}


void SceneCloth::mousePressed(const QMouseEvent* e, const Camera& cam)
{
    grabX = e->pos().x();
    grabY = e->pos().y();

    if (!(e->modifiers() & Qt::ControlModifier)) {

        Vec3 rayDir    = cam.getRayDir(grabX, grabY);
        Vec3 rayOrigin = cam.getPos();

        selectedParticle = -1;
        for (int i = 0; i < numParticles; i++) {
            Particle* p =  system.getParticle(i);
            bool hitParticle = checkIntersectionWithParticle(rayOrigin, rayDir, p);

            if (hitParticle)
                selectedParticle = i;
        }

        if (selectedParticle >= 0) {
            cursorWorldPos = system.getParticle(selectedParticle)->pos;
        }
    }
}


void SceneCloth::mouseMoved(const QMouseEvent* e, const Camera& cam)
{
    int dx = e->pos().x() - grabX;
    int dy = e->pos().y() - grabY;
    grabX = e->pos().x();
    grabY = e->pos().y();

    Vec3 disp = cam.worldSpaceDisplacement(dx, -dy, cam.getEyeDistance());

    if (e->modifiers() & Qt::ControlModifier) {
        // Move props
        Vec3 rayDir = cam.getRayDir(grabX, grabY);
        Vec3 rayOrigin = cam.getPos();

        bool hitAABB   = checkIntersectionWithAABB(rayOrigin, rayDir, colliderPropBox);
        bool hitSphere = checkIntersectionWithSphere(rayOrigin, rayDir, colliderPropSphere);

        if (hitAABB)
            colliderPropBox.Move(disp);
        if (hitSphere)
            colliderPropSphere.Move(disp);
    }
    else if (e->modifiers() & Qt::ShiftModifier) {

    }
    else {
        if (selectedParticle >= 0) {
            double d = -(system.getParticle(selectedParticle)->pos - cam.getPos()).dot(cam.zAxis());
            Vec3 disp = cam.worldSpaceDisplacement(dx, -dy, d);
            cursorWorldPos += disp;
        }
    }
}



void SceneCloth::mouseReleased(const QMouseEvent*, const Camera&)
{
    selectedParticle = -1;
}

void SceneCloth::keyPressed(const QKeyEvent* e, const Camera&)
{
    if (selectedParticle >= 0 && e->key() == Qt::Key_F) {
        fixedParticle[selectedParticle] = true;
        Particle* p = system.getParticle(selectedParticle);
        p->prevPos = p->pos;
        p->vel = Vec3(0,0,0);
        p->force = Vec3(0,0,0);
    }
}


// Mouse raycast interaction
bool SceneCloth::checkIntersectionWithParticle(const Vec3& rayOrigin, const Vec3& rayDir, Particle* p)
{
    Vec3 radius = Vec3(p->radius, p->radius, p->radius);
    Eigen::Array3d tMinArray = (p->pos - radius - rayOrigin).array() / rayDir.array();
    Eigen::Array3d tMaxArray = (p->pos + radius - rayOrigin).array() / rayDir.array();

    Eigen::Array3d t1Array = (tMinArray < tMaxArray).select(tMinArray, tMaxArray);  // min of the 2
    Eigen::Array3d t2Array = (tMinArray < tMaxArray).select(tMaxArray, tMinArray);  // max of the 2


    float tNear = std::max(std::max(t1Array.x(), t1Array.y()), t1Array.z());    // max of the 3
    float tFar =  std::min(std::min(t2Array.x(), t2Array.y()), t2Array.z());    // min of the 3

    if (tNear <= tFar)
        return true;

    return false;
}

bool SceneCloth::checkIntersectionWithAABB(const Vec3& rayOrigin, const Vec3& rayDir, const ColliderAABB& collider)
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

bool SceneCloth::checkIntersectionWithSphere(const Vec3& rayOrigin, const Vec3& rayDir, const ColliderSphere& collider)
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



// Compute Cloth
void SceneCloth::computeCloth(int particles_width, int particles_height,
                        QVector<Particle*>& sceneParticles, QVector<ForceSpring*>& sceneForceSprings,
                        double clothWidth, double clothHeight, double particleRadius) {

    // Grid edges
    double edgeX = clothWidth / particles_width;
    double edgeY = clothHeight / particles_height;

    // Create particles in a grid
    for (int i = 0; i < particles_width; i++) {
        for (int j = 0; j < particles_height; j++) {
            double tx = i * edgeX - 0.5 * clothWidth;
            double ty = j * edgeY - 0.5 * clothHeight;
            Vec3 pos = Vec3(ty + edgeY, 70 - edgeX - tx, 0);

            Particle* newParticle = new Particle();
            newParticle->id = i * particles_height + j;
            newParticle->pos = pos;
            newParticle->prevPos = pos;
            newParticle->vel = Vec3(0, 0, 0);
            newParticle->mass = 1;
            newParticle->radius = particleRadius;
            newParticle->color = Vec3(235 / 255.0, 51 / 255.0, 36 / 255.0);
            newParticle->fixed = false;

            // Locking the first and last particle in the first row (optional)
            if ((i == 0 && j == 0) || (i == 0 && j == particles_height - 1)) {
                newParticle->fixed = true;
            }

            sceneParticles.push_back(newParticle);
        }
    }

    // Create springs for stretch, shear, and bend
    for (int i=0; i<particles_width; i++) {
        for (int j=0; j<particles_height; j++) {

            int idx = i * particles_height + j;
            Particle* p0 = sceneParticles[idx]; // current particle

            // Stretch
            // ------------

            // Right
            if (i < particles_width - 1) {
                int rightNeighborIdx = (i + 1) * particles_height + j;
                Particle* p1 = sceneParticles[rightNeighborIdx];
                ForceSpring* spring = new ForceSpring(p0, p1);
                spring->addInfluencedParticle(p0);
                spring->addInfluencedParticle(p1);
                sceneForceSprings.push_back(spring);
            }
            // Down
            if (j < particles_height - 1) {
                int downNeighborIdx = i * particles_height + (j + 1);
                Particle* p1 = sceneParticles[downNeighborIdx];
                ForceSpring* spring = new ForceSpring(p0, p1);
                spring->addInfluencedParticle(p0);
                spring->addInfluencedParticle(p1);
                sceneForceSprings.push_back(spring);
            }

            // Shear
            // ------------

            // Down-Right
            if (i < particles_width - 1 && j < particles_height - 1) {
                int downRightNeighborIdx = (i + 1) * particles_height + (j + 1);
                Particle* p1 = sceneParticles[downRightNeighborIdx];
                ForceSpring* spring = new ForceSpring(p0, p1);
                spring->addInfluencedParticle(p0);
                spring->addInfluencedParticle(p1);
                sceneForceSprings.push_back(spring);
            }

            // Down-Left
            if (i > 0 && j < particles_height - 1) {
                int downLeftNeighborIdx = (i - 1) * particles_height + (j + 1);
                Particle* p1 = sceneParticles[downLeftNeighborIdx];
                ForceSpring* spring = new ForceSpring(p0, p1);
                spring->addInfluencedParticle(p0);
                spring->addInfluencedParticle(p1);
                sceneForceSprings.push_back(spring);
            }

            // Bend
            // ------------

            // Right (2 cells away)
            if (i < particles_width - 2) {
                int rightNeighborIdx = (i + 2) * particles_height + j;
                Particle* p1 = sceneParticles[rightNeighborIdx];
                ForceSpring* spring = new ForceSpring(p0, p1);
                spring->addInfluencedParticle(p0);
                spring->addInfluencedParticle(p1);
                sceneForceSprings.push_back(spring);
            }

            // Down (2 cells away)
            if (j < particles_height - 2) {
                int downNeighborIdx = i * particles_height + (j + 2);
                Particle* p1 = sceneParticles[downNeighborIdx];
                ForceSpring* spring = new ForceSpring(p0, p1);
                spring->addInfluencedParticle(p0);
                spring->addInfluencedParticle(p1);
                sceneForceSprings.push_back(spring);
            }
        }
    }
}

