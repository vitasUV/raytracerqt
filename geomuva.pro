TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

include(deployment.pri)
qtcAddDeployment()

HEADERS += \
    geometry.h

DISTFILES +=

