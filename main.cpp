#include "MainWindow.h"
#include "Model.h"
#include "Visualizer.h"
#include "PowerSearch.h"

#include <QApplication>
#include <QtGui>
#include <QThread>
#include <QObject>



int main(int argc, char *argv[])
{
    PowerSearch ps(argc, argv);



    return ps.exec();
}
