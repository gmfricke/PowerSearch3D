#include "Worker.h"
#include "Model.h"
#include <iostream>

using namespace std;

Worker::Worker(QObject *parent) :
QThread(parent)
{
    pause = false;
    generate_paths = true;
    restart = false;
    _abort = false;
}


Worker::~Worker()
{
}

void Worker::abort()
{
mutex.lock();
_abort = true;
condition.wakeOne();
mutex.unlock();

wait();
}

void Worker::startWork()
{
QMutexLocker locker(&mutex);

if (!isRunning()) {
start(LowPriority);
} else {
restart = true;
condition.wakeOne();
}
}

void Worker::setModel(Model* m)
{
    model = m;
}

void Worker::run()
{

//forever {
    // Perform central task

    cout << "Worker: running in thread " <<  pthread_self() << endl;

    /*
    pause = true;
    while (pause)
    {
        cout << "Worker paused" << endl;
        usleep(100);
    }
    cout << "Worker woke up" << endl;
*/
  model->doAnalysis();

//if(!restart) emit workFinished();

  //cout << "Worker: work finished" << endl;
  //emit workFinished();
  //quit();
 // cout << "Worker: quit()" << endl;

// Abort thread if requested
//if(_abort) return;

// Restart when new central task is given
//mutex.lock();
//if(!restart) condition.wait(&mutex);
//restart = false;
//mutex.unlock();
//}
}

void Worker::setGeneratePaths( bool gp )
{
    generate_paths = gp;
}

void Worker::wakeup()
{
    cout << "Worker: wakeup called" << endl;
    pause = false;
}
