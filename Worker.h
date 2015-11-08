#ifndef WORKER_H
#define WORKER_H

#include "Model.h"
#include <QThread>
#include <QMutex>
#include <QWaitCondition>

class Model;

class Worker : public QThread
{
Q_OBJECT

public:
explicit Worker(QObject *parent = 0);
    void abort();
    void setModel(Model* m);
    void setGeneratePaths(bool gp);
    ~Worker();

public slots:
    void startWork();
    void wakeup();

signals:
void reportProgress(int pProg, QString message);
void workFinished();

protected:
void run();


QMutex mutex;
QWaitCondition condition;
bool restart;
bool _abort;
bool generate_paths;

private:

bool pause;
Model* model;
};

#endif // WORKER_H
