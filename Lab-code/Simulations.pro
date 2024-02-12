QT       += core gui opengl
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++11

INCLUDEPATH += code
INCLUDEPATH += extlibs

VPATH += code

SOURCES += \
    code/camera.cpp \
    code/colliders.cpp \
    code/forces.cpp \
    code/glutils.cpp \
    code/glwidget.cpp \
    code/integrators.cpp \
    code/main.cpp \
    code/mainwindow.cpp \
    code/model.cpp \
    code/particlesystem.cpp \
    code/scenefluid.cpp \
    code/scenefountain.cpp \
    code/sceneprojectiles.cpp \
    code/scenecloth.cpp \
    code/widgetfluid.cpp \
    code/widgetfountain.cpp \
    code/widgetprojectiles.cpp \
    code/widgetcloth.cpp \

HEADERS += \
    code/camera.h \
    code/colliders.h \
    code/defines.h \
    code/forces.h \
    code/glutils.h \
    code/glwidget.h \
    code/integrators.h \
    code/mainwindow.h \
    code/model.h \
    code/particle.h \
    code/particlesystem.h \
    code/scene.h \
    code/scenefluid.h \
    code/scenefountain.h \
    code/sceneprojectiles.h \
    code/scenecloth.h \
    code/spatialhash.h \
    code/widgetfluid.h \
    code/widgetfountain.h \
    code/widgetprojectiles.h \
    code/widgetcloth.h \

FORMS += \
    forms/mainwindow.ui \
    forms/widgetfluid.ui \
    forms/widgetfountain.ui \
    forms/widgetprojectiles.ui \
    forms/widgetcloth.ui \

RESOURCES += shaders.qrc
