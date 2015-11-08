#-------------------------------------------------
#
# Project created by QtCreator 2013-03-30T12:17:29
#
#-------------------------------------------------

QT       += core gui opengl printsupport
QMAKE_LFLAGS += -v

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = PowerSearch3Dv8
TEMPLATE = app


INCLUDEPATH += /Users/matthew/Dropbox/ToR/TCellMotility/Code/Simulation/PowerSearch3Dv8/ALGLIB/
SOURCES += /Users/matthew/Dropbox/ToR/TCellMotility/Code/Simulation/PowerSearch3Dv8/ALGLIB/*.cpp

SOURCES += main.cpp\
    Model.cpp \
    SearchSpace.cpp \
    Visualizer.cpp \
    Agent.cpp \
    Searcher.cpp \
    Target.cpp \
    Coordinate.cpp \
    rand_gen.cpp \
    MainWindow.cpp \
    RectangularPrism.cpp \
    Sphere.cpp \
    Worker.cpp \
    PowerSearch.cpp\
    Histogram.cpp\
    QCustomPlot.cpp


HEADERS  += \
    Model.h \
    SearchSpace.h \
    Visualizer.h \
    Agent.h \
    Searcher.h \
    Target.h \
    Coordinate.h \
    rand_gen.h \
    MainWindow.h \
    RectangularPrism.h \
    Sphere.h \
    Worker.h \
    PowerSearch.h \
    Histogram.h\
    QCustomPlot.h \
    QDebugStream.h


FORMS    += mainwindow.ui

mac: LIBS += -framework GLUT
INCLUDEPATH += /opt/local/include/
DEPENDPATH += /opt/local/include/

ICON = powersearch.icns
QMAKE_INFO_PLIST = powersearch.plist # Use custom .plist to take advantage of Retina display

#unix: LIBS += -L/usr/lib/ -lboost_random-mt

INCLUDEPATH += /usr/include
DEPENDPATH += /usr/include

#unix: PRE_TARGETDEPS += /usr/lib/libboost_random-mt.a

unix:!mac: LIBS += -lglut -lGLU

INCLUDEPATH += -L/usr/lib/
