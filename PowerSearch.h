#ifndef POWERSEARCH_H
#define POWERSEARCH_H

#include "Model.h"
#include "Visualizer.h"
#include "MainWindow.h"

#include <QtGui>
#include <QThread>

#include <QApplication>

class PowerSearch : public QApplication
{
    Q_OBJECT
public:
    explicit PowerSearch(int argc, char *argv[]);

signals:

public slots:
    void quit();
private:
    Model* model;
    Visualizer* vis;
    MainWindow* window;

};

#endif // POWERSEARCH_H
