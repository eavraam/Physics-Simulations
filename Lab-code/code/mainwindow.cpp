#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "sceneprojectiles.h"
#include "scenefountain.h"
#include "scenecloth.h"
#include "scenefluid.h"
//
//
//
//
//

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    connect(ui->btnSimstep, SIGNAL(clicked()), ui->openGLWidget, SLOT(doSimStep()));
    connect(ui->btnSimloop, SIGNAL(clicked()), ui->openGLWidget, SLOT(doSimLoop()));
    connect(ui->btnPause,   SIGNAL(clicked()), ui->openGLWidget, SLOT(pauseSim()));
    connect(ui->btnReset,   SIGNAL(clicked()), ui->openGLWidget, SLOT(resetSim()));
    connect(ui->timestep,   SIGNAL(valueChanged(double)), ui->openGLWidget, SLOT(setTimeStep(double)));
    ui->openGLWidget->setTimeStep(ui->timestep->value());

    connect(ui->blackHoleToggle,   SIGNAL(stateChanged(int)), ui->openGLWidget, SLOT(setBlackHole(int)));

    connect(ui->kElastic,   SIGNAL(valueChanged(double)), ui->openGLWidget, SLOT(setElasticity(double)));
    connect(ui->kFriction,   SIGNAL(valueChanged(double)), ui->openGLWidget, SLOT(setFriction(double)));

    connect(ui->drag_none,   SIGNAL(clicked()), ui->openGLWidget, SLOT(setNoDrag()));
    connect(ui->drag_linear,   SIGNAL(clicked()), ui->openGLWidget, SLOT(setLinearDrag()));
    connect(ui->drag_quadratic,   SIGNAL(clicked()), ui->openGLWidget, SLOT(setQuadraticDrag()));

    connect(ui->actionCameraReset, SIGNAL(triggered()), ui->openGLWidget, SLOT(resetCamera()));
    connect(ui->actionCameraX, SIGNAL(triggered()), ui->openGLWidget, SLOT(cameraViewX()));
    connect(ui->actionCameraY, SIGNAL(triggered()), ui->openGLWidget, SLOT(cameraViewY()));
    connect(ui->actionCameraZ, SIGNAL(triggered()), ui->openGLWidget, SLOT(cameraViewZ()));

    connect(ui->actionSceneFountain,   &QAction::triggered, this, [=](void){ changeScene(new SceneFountain()); });
    connect(ui->actionSceneProjectiles,&QAction::triggered, this, [=](void){ changeScene(new SceneProjectiles()); });
    connect(ui->actionSceneCloth,      &QAction::triggered, this, [=](void){ changeScene(new SceneCloth()); });
    connect(ui->actionSceneFluid,      &QAction::triggered, this, [=](void){ changeScene(new SceneFluid()); });
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::changeScene(Scene* sc)
{
    // First set scene to glwidget. This will trigger deleting old scene and its associated widget.
    // Otherwise, qDeleteAll deletes the widget first, and then the scene destructor crashes
    ui->openGLWidget->setScene(sc);
    qDeleteAll(ui->controlsScene->findChildren<QWidget *>(QString(), Qt::FindDirectChildrenOnly));
    ui->controlsScene->layout()->addWidget(sc->sceneUI());
}
