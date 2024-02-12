#include "sceneprojectiles.h"
#include "glutils.h"
#include "model.h"
#include <QMatrix4x4>
#include <iostream>


SceneProjectiles::SceneProjectiles() {
    widget = new WidgetProjectiles();

    const std::vector<std::string> solvers = {
        "Euler", "Symplectic Euler", "Verlet", "Midpoint",
        "RK-2", "Velocity Verlet", "RK-4"
    };
    widget->setSolverTypes(solvers);
    widget->setSolver1(0); // Euler
    widget->setSolver2(2); // Symplectic Euler
}

SceneProjectiles::~SceneProjectiles() {
    if (widget)      delete widget;
    if (shaderPhong) delete shaderPhong;
    if (shaderLines) delete shaderLines;
    if (vaoSphere)   delete vaoSphere;
    if (vaoFloor)    delete vaoFloor;
    if (vaoCube)     delete vaoCube;
    if (vaoTrajectory) delete vaoTrajectory;
    if (vboTrajectoryPoints) delete vboTrajectoryPoints;
    if (integrator1) delete integrator1;
    if (integrator2) delete integrator2;
    systemAnalytic.deleteParticles();
    systemAnalytic.deleteForces();
    systemNumerical1.deleteParticles();
    systemNumerical1.deleteForces();
    systemNumerical2.deleteParticles();
    systemNumerical2.deleteForces();
}

void SceneProjectiles::initialize(double dt, int dr, double kEl, double kFr, int bht) {

    // Initialize drag, elasticity, friction and black hole
    timeStep = dt;
    drag = dr;
    kElastic = kEl;
    kFriction = kFr;
    blackHoleToggle = bht;


    std::cout << "Initial kElastic: " << kElastic << std::endl;
    std::cout << "Initial kFriction: " << kFriction << std::endl;

    /*
     * RENDER INITS
     */

    // load shader programs
    shaderPhong = glutils::loadShaderProgram(":/shaders/phong.vert", ":/shaders/phong.frag");
    shaderLines = glutils::loadShaderProgram(":/shaders/lines.vert", ":/shaders/lines.geom", ":/shaders/lines.frag");

    // create sphere VAO
    Model sphere = Model::createIcosphere(3);
    vaoSphere = glutils::createVAO(shaderPhong, &sphere);
    numSphereFaces = sphere.numFaces();
    glutils::checkGLError();

    // create floor VAO
    Model quad = Model::createQuad();
    vaoFloor = glutils::createVAO(shaderPhong, &quad);
    glutils::checkGLError();

    // create cube VAO
    Model cube = Model::createCube();
    vaoCube = glutils::createVAO(shaderPhong, &cube);
    glutils::checkGLError();

    // initialize the buffer we'll use to render trajectories
    vaoTrajectory = new QOpenGLVertexArrayObject();
    vaoTrajectory->create();
    vaoTrajectory->bind();
    vboTrajectoryPoints = new QOpenGLBuffer(QOpenGLBuffer::Type::VertexBuffer);
    vboTrajectoryPoints->create();
    vboTrajectoryPoints->bind();
    vboTrajectoryPoints->setUsagePattern(QOpenGLBuffer::UsagePattern::DynamicDraw);
    vboTrajectoryPoints->allocate(MAX_TRAJ_POINTS*3*sizeof(float));
    shaderLines->setAttributeBuffer("vertex", GL_FLOAT, 0, 3, 0);
    shaderLines->enableAttributeArray("vertex");
    vaoTrajectory->release();


    /*
     * PHYSICS INITS
     */

    // create the different systems, with one particle each
    systemAnalytic.addParticle(new Particle(Vec3( 0, 0, 0), Vec3(0,0,0), 1));
    systemAnalytic.getParticle(0)->color = Vec3(0, 0.5, 0);
    systemAnalytic.getParticle(0)->radius = 2;
    systemNumerical1.addParticle(new Particle(Vec3( 0, 0, 15), Vec3(0,0,0), 1));
    systemNumerical1.getParticle(0)->color = Vec3(0.5, 0, 0);
    systemNumerical1.getParticle(0)->radius = 2;
    systemNumerical2.addParticle(new Particle(Vec3( 0, 0, -15), Vec3(0,0,0), 1));
    systemNumerical2.getParticle(0)->color = Vec3(0, 0, 0.5);
    systemNumerical2.getParticle(0)->radius = 2;

    // prevpos
    // I set these in the reset function too, and glwidget>setScene calls both initialize() and reset() functions.
    // Hence, these 3 lines are not necessary, but I leave them for readability and in case of any glwidget change.
    systemAnalytic.getParticle(0)->prevPos = systemAnalytic.getParticle(0)->pos - timeStep*systemAnalytic.getParticle(0)->vel;
    systemNumerical1.getParticle(0)->prevPos = systemNumerical1.getParticle(0)->pos - timeStep*systemNumerical1.getParticle(0)->vel;
    systemNumerical2.getParticle(0)->prevPos = systemNumerical2.getParticle(0)->pos - timeStep*systemNumerical2.getParticle(0)->vel;

    // floor
    colliderFloor.setPlane(Vec3(0, 1, 0), 0);


    // Forces
    // we need to create one per system to assign its particle
    // reset them in the init, just in case
    systemNumerical1.clearForces();
    systemNumerical2.clearForces();

    // gravity - ALWAYS
    fGravity1 = new ForceConstAcceleration(Vec3(0, -gravityAccel, 0));
    fGravity1->addInfluencedParticle(systemNumerical1.getParticle(0));
    systemNumerical1.addForce(fGravity1);

    fGravity2 = new ForceConstAcceleration(Vec3(0, -gravityAccel, 0));
    fGravity2->addInfluencedParticle(systemNumerical2.getParticle(0));
    systemNumerical2.addForce(fGravity2);

    // drag is set to 0 by default, but we need these initializations
    // in case the user manually changes the default value
    // drag linear
    if (drag == 1)
    {
        fDragLinear1 = new ForceDragLinear();
        fDragLinear1->addInfluencedParticle(systemNumerical1.getParticle(0));
        systemNumerical1.addForce(fDragLinear1);

        fDragLinear2 = new ForceDragLinear();
        fDragLinear2->addInfluencedParticle(systemNumerical2.getParticle(0));
        systemNumerical2.addForce(fDragLinear2);
    }
    // drag quadratic
    else if (drag == 2)
    {
        fDragQuadratic1 = new ForceDragQuadratic(Vec3());
        fDragQuadratic1->addInfluencedParticle(systemNumerical1.getParticle(0));
        systemNumerical1.addForce(fDragQuadratic1);

        fDragQuadratic2 = new ForceDragQuadratic(Vec3());
        fDragQuadratic2->addInfluencedParticle(systemNumerical2.getParticle(0));
        systemNumerical2.addForce(fDragQuadratic2);
    }
}


Integrator* createIntegrator(int type) {
    switch(type) {
        case 0: return new IntegratorEuler();
        case 1: return new IntegratorSymplecticEuler();
        case 2: return new IntegratorVerlet();
        case 3: return new IntegratorMidpoint();
        case 4: return new IntegratorRK2();
        case 5: return new IntegratorVelocityVerlet();
        case 6: return new IntegratorRK4();
        default: return nullptr;
    }
}

void SceneProjectiles::reset(double dt, int dr, double kEl, double kFr, int bht) {

    // Reset drag, elasticity, friction and black hole
    timeStep = dt;
    drag = dr;
    kElastic = kEl;
    kFriction = kFr;
    blackHoleToggle = bht;

    std::cout << "Reset kElastic to: " << kElastic << std::endl;
    std::cout << "Reset kFriction to: " << kFriction << std::endl;
    std::cout << "Drag reset to: " << drag << std::endl;
    std::cout << "--------------- \n" << drag << std::endl;

    // get new interface values
    shotHeight   = widget->getHeight();
    shotAngle    = Math::toRad(widget->getAngle());
    shotSpeed    = widget->getSpeed();
    gravityAccel = widget->getGravity();

    // integrators
    if (integrator1) delete integrator1;
    if (integrator2) delete integrator2;
    integrator1 = createIntegrator(widget->getSolver1());
    integrator2 = createIntegrator(widget->getSolver2());

    // update initial particle positions
    const double zdist = 15;
    systemAnalytic.getParticle(0)->pos = Vec3(0, shotHeight, 0);
    systemAnalytic.getParticle(0)->vel = shotSpeed*Vec3(std::cos(shotAngle), std::sin(shotAngle), 0);
    systemAnalytic.getParticle(0)->prevPos = systemAnalytic.getParticle(0)->pos - timeStep*systemAnalytic.getParticle(0)->vel;
    systemNumerical1.getParticle(0)->pos = Vec3(0, shotHeight, zdist);
    systemNumerical1.getParticle(0)->vel = shotSpeed*Vec3(std::cos(shotAngle), std::sin(shotAngle), 0);
    systemNumerical1.getParticle(0)->prevPos = systemNumerical1.getParticle(0)->pos - timeStep*systemNumerical1.getParticle(0)->vel;
    systemNumerical2.getParticle(0)->pos = Vec3(0, shotHeight, -zdist);
    systemNumerical2.getParticle(0)->vel = shotSpeed*Vec3(std::cos(shotAngle), std::sin(shotAngle), 0);
    systemNumerical2.getParticle(0)->prevPos = systemNumerical2.getParticle(0)->pos - timeStep*systemNumerical2.getParticle(0)->vel;


    // Forces
    // reset system forces
    systemNumerical1.clearForces();
    systemNumerical2.clearForces();

    // update gravity accelerations
    // since we clear the forces, we need to reinitialize, or else the FIRST gravity will work correctly
    // then switching to linear-quadratic will work, and returning to only-gravity will have
    // persistent values from the others.
    fGravity1 = new ForceConstAcceleration(Vec3(0, -gravityAccel, 0));
    fGravity1->addInfluencedParticle(systemNumerical1.getParticle(0));
    systemNumerical1.addForce(fGravity1);

    fGravity2 = new ForceConstAcceleration(Vec3(0, -gravityAccel, 0));
    fGravity2->addInfluencedParticle(systemNumerical2.getParticle(0));
    systemNumerical2.addForce(fGravity2);

    // drag is set to 0 by default, so we need these initializations
    // in case the user restarts with drag = 1 or 2, since they were not initialized
    // update drag linear accelerations
    if (drag == 1)
    {
        fDragLinear1 = new ForceDragLinear();
        fDragLinear1->addInfluencedParticle(systemNumerical1.getParticle(0));
        systemNumerical1.addForce(fDragLinear1);

        fDragLinear2 = new ForceDragLinear();
        fDragLinear2->addInfluencedParticle(systemNumerical2.getParticle(0));
        systemNumerical2.addForce(fDragLinear2);
    }
    // drag quadratic
    else if (drag == 2)
    {
        fDragQuadratic1 = new ForceDragQuadratic(Vec3());
        fDragQuadratic1->addInfluencedParticle(systemNumerical1.getParticle(0));
        systemNumerical1.addForce(fDragQuadratic1);

        fDragQuadratic2 = new ForceDragQuadratic(Vec3());
        fDragQuadratic2->addInfluencedParticle(systemNumerical2.getParticle(0));
        systemNumerical2.addForce(fDragQuadratic2);
    }

    // update system forces
    systemNumerical1.updateForces();
    systemNumerical2.updateForces();

    // trajectories
    trajectoryAnalytic.clear();
    trajectoryAnalytic.push_back(systemAnalytic.getParticle(0)->pos);
    trajectoryNumerical1.clear();
    trajectoryNumerical1.push_back(systemNumerical1.getParticle(0)->pos);
    trajectoryNumerical2.clear();
    trajectoryNumerical2.push_back(systemNumerical2.getParticle(0)->pos);

    // put particles to run
    system1active = true;
    system2active = true;

    // reset timer
    time = 0;
}


void SceneProjectiles::update(double dt) {

    // total ellapsed time (needed for analytic solution)
    time += dt;

    // ANALYTIC: projectile motion equations until we reach the ground
    Particle* p = systemAnalytic.getParticle(0);
    double vy0 = shotSpeed*std::sin(shotAngle);
    double tGround = (vy0 + std::sqrt(vy0*vy0 + 2*gravityAccel*shotHeight))/gravityAccel;
    if (time - dt <= tGround) {
        double t = std::min(time, tGround);
        p->pos[0] = t * shotSpeed * std::cos(shotAngle);
        p->pos[1] = shotHeight + t*vy0 - 0.5*gravityAccel*t*t;
        p->vel    = Vec3(shotSpeed*std::cos(shotAngle),
                         shotSpeed*std::sin(shotAngle) - gravityAccel*t, 0);

        trajectoryAnalytic.push_back(p->pos);
        if (trajectoryAnalytic.size() > MAX_TRAJ_POINTS) trajectoryAnalytic.pop_front();
    }


    // NUMERICAL INTEGRATORS:

    if (system1active) {
        // integration step
        integrator1->step(systemNumerical1, dt);

        // collision test
        Particle* p = systemNumerical1.getParticle(0);

        if (colliderFloor.testCollision(p, dt)) {
            // std::cout << "Particle 1: passed the collision test..." << std::endl;
            colliderFloor.resolveCollision(p, kElastic, kFriction);

            // Update prevPos for the Verlet solver
            if (widget->getSolver1() == 2) {
                p->prevPos = p->pos - (p->vel * dt);
            }
        }

        // record trajectory
        trajectoryNumerical1.push_back(p->pos);
        if (trajectoryNumerical1.size() > MAX_TRAJ_POINTS) {
            trajectoryNumerical1.pop_front();
        }

    }

    if (system2active) {
        // integration step
        integrator2->step(systemNumerical2, dt);

        // collision test
        Particle* p = systemNumerical2.getParticle(0);

        if (colliderFloor.testCollision(p, dt)) {
            // std::cout << "Particle 2: passed the collision test..." << std::endl;
            colliderFloor.resolveCollision(p, kElastic, kFriction);

            // Update prevPos for the Verlet solver
            if (widget->getSolver2() == 2) {
                p->prevPos = p->pos - (p->vel * dt);
            }
        }

        // record trajectory
        trajectoryNumerical2.push_back(p->pos);
        if (trajectoryNumerical2.size() > MAX_TRAJ_POINTS) {
            trajectoryNumerical2.pop_front();
        }
    }


    // print particle heights
    // std::cout << "time: " << time << std::endl;
    // std::cout << "analytic sol: " << systemAnalytic.getParticle(0)->pos[1] << std::endl;
    // std::cout << "numerical 1: " << systemNumerical1.getParticle(0)->pos[1] << std::endl;
    // std::cout << "numerical 2: " << systemNumerical2.getParticle(0)->pos[1] << std::endl;
    // std::cout << std::endl;
}


void SceneProjectiles::paint(const Camera& camera) {

    // pointer to current context OpenGL functions
    QOpenGLFunctions *glFuncs = QOpenGLContext::currentContext()->functions();

    // start using phong shader
    shaderPhong->bind();

    // camera matrices
    QMatrix4x4 camProj = camera.getPerspectiveMatrix();
    QMatrix4x4 camView = camera.getViewMatrix();
    shaderPhong->setUniformValue("ProjMatrix", camProj);
    shaderPhong->setUniformValue("ViewMatrix", camView);

    // lighting
    const int numLights = 2;
    const QVector3D lightPosWorld[numLights] = {QVector3D(1000,1000,1000), QVector3D(-1000,1000,-1000)};
    const QVector3D lightColor[numLights] = {QVector3D(1,1,1), QVector3D(1,1,1)};
    QVector3D lightPosCam[numLights];
    for (int i = 0; i < numLights; i++) {
        lightPosCam[i] = camView * lightPosWorld[i];
    }
    shaderPhong->setUniformValue("numLights", numLights);
    shaderPhong->setUniformValueArray("lightPos", lightPosCam, numLights);
    shaderPhong->setUniformValueArray("lightColor", lightColor, numLights);

    // draw floor
    vaoFloor->bind();
    QMatrix4x4 modelMat;
    modelMat.translate(150, 0, 0);
    modelMat.scale(200, 1, 50);
    shaderPhong->setUniformValue("ModelMatrix", modelMat);
    shaderPhong->setUniformValue("matdiff", 0.7f, 0.7f, 0.7f);
    shaderPhong->setUniformValue("matspec", 0.0f, 0.0f, 0.0f);
    shaderPhong->setUniformValue("matshin", 0.0f);
    glFuncs->glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

    // draw cube
    vaoCube->bind();
    modelMat = QMatrix4x4();
    modelMat.translate(-1, 0.5*shotHeight - 1, 0);
    modelMat.scale(2, 0.5*shotHeight - 1, 25.0f);
    shaderPhong->setUniformValue("ModelMatrix", modelMat);
    glFuncs->glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);

    // draw the different spheres
    vaoSphere->bind();
    const Particle* particles[3] = { systemAnalytic.getParticle(0),
                                     systemNumerical1.getParticle(0),
                                     systemNumerical2.getParticle(0) };
    for (const Particle* particle : particles) {
        Vec3   p = particle->pos;
        Vec3   c = particle->color;
        double r = particle->radius;

        modelMat = QMatrix4x4();
        modelMat.translate(p[0], p[1], widget->renderSameZ() ? 0 : p[2]);
        modelMat.scale(r);
        shaderPhong->setUniformValue("ModelMatrix", modelMat);

        shaderPhong->setUniformValue("matdiff", GLfloat(c[0]), GLfloat(c[1]), GLfloat(c[2]));
        shaderPhong->setUniformValue("matspec", 1.0f, 1.0f, 1.0f);
        shaderPhong->setUniformValue("matshin", 100.f);

        glFuncs->glDrawElements(GL_TRIANGLES, 3*numSphereFaces, GL_UNSIGNED_INT, 0);
    }
    vaoSphere->release();

    // we're done with this shader
    shaderPhong->release();

    // do we need to draw trajectories?
    if (widget->renderTrajectory()) {

        shaderLines->bind();
        shaderLines->setUniformValue("ProjMatrix", camProj);
        shaderLines->setUniformValue("ViewMatrix", camView);
        shaderLines->setUniformValue("radius", float(0.25));
        shaderLines->setUniformValue("shading", false);

        updateTrajectoryCoordsBuffer(trajectoryAnalytic, widget->renderSameZ());
        Vec3 c = systemAnalytic.getParticle(0)->color;
        shaderLines->setUniformValue("matdiff", GLfloat(c[0]), GLfloat(c[1]), GLfloat(c[2]));
        vaoTrajectory->bind();
        glFuncs->glDrawArrays(GL_LINE_STRIP, 0, std::min(static_cast<unsigned int>(trajectoryAnalytic.size()),
                                                         MAX_TRAJ_POINTS));
        vaoTrajectory->release();

        updateTrajectoryCoordsBuffer(trajectoryNumerical1, widget->renderSameZ());
        c = systemNumerical1.getParticle(0)->color;
        shaderLines->setUniformValue("matdiff", GLfloat(c[0]), GLfloat(c[1]), GLfloat(c[2]));
        vaoTrajectory->bind();
        glFuncs->glDrawArrays(GL_LINE_STRIP, 0, std::min(static_cast<unsigned int>(trajectoryNumerical1.size()),
                                                         MAX_TRAJ_POINTS));
        vaoTrajectory->release();

        updateTrajectoryCoordsBuffer(trajectoryNumerical2, widget->renderSameZ());
        c = systemNumerical2.getParticle(0)->color;
        shaderLines->setUniformValue("matdiff", GLfloat(c[0]), GLfloat(c[1]), GLfloat(c[2]));
        vaoTrajectory->bind();
        glFuncs->glDrawArrays(GL_LINE_STRIP, 0, std::min(static_cast<unsigned int>(trajectoryNumerical2.size()),
                                                         MAX_TRAJ_POINTS));
        vaoTrajectory->release();

        shaderLines->release();
    }
}


void SceneProjectiles::updateTrajectoryCoordsBuffer(const std::list<Vec3>& trajectory, bool sameZ) {
    vboTrajectoryPoints->bind();
    float* pos = new float[3*trajectory.size()];
    unsigned int i = 0;
    for (auto it = trajectory.begin(); it != trajectory.end(); it++) {
        pos[3*i  ] = it->x();
        pos[3*i+1] = it->y();
        pos[3*i+2] = sameZ ? 0 : it->z();
        i++;
    }
    void* bufptr = vboTrajectoryPoints->mapRange(0, 3*trajectory.size()*sizeof(float),
                                                 QOpenGLBuffer::RangeInvalidateBuffer | QOpenGLBuffer::RangeWrite);
    memcpy(bufptr, (void*)(pos), 3*trajectory.size()*sizeof(float));
    vboTrajectoryPoints->unmap();
    vboTrajectoryPoints->release();
    delete[] pos;
}
