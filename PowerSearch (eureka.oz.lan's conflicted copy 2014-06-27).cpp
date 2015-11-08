#include "PowerSearch.h"
#include "iostream"
#include "QDebugStream.h"

using namespace std;

PowerSearch::PowerSearch(int argc, char *argv[]) : QApplication(argc, argv)
{
    // Call quit when the last window is closed
    QObject::connect( this, SIGNAL(lastWindowClosed()), this, SLOT(quit()) );

    qsrand(QTime(0,0,0).secsTo(QTime::currentTime()));

    window = new MainWindow();

    float spr= 5.00, cr = 40.00;
    int nc = 20, ntc = 10, sp = 4;
    bool t = true;

    string app_path = applicationDirPath().toStdString();
    model = new Model(spr, cr, nc, ntc, sp, t, app_path);
    cout << "Determined application path to be " << app_path << endl;
    vis = window->getVisualizer();
    vis->setModel(model);

    model->window = window;

    //QThread *visualizerThread = new QThread(vis);
    qRegisterMetaType<MyArray>("MyArray");
    qRegisterMetaType<FloatVector>("FloatVector");

    //QObject::connect(model, SIGNAL(updateNeeded()), window, SLOT(update()));
    QObject::connect(model, SIGNAL(updateNeeded()), vis, SLOT(updateGL()));
    QObject::connect(model, SIGNAL(parseError(QString)), vis, SLOT(displayParseError(QString)));
    QObject::connect(model, SIGNAL(updateTimeEfficiencyMean(float)), window, SLOT(setTimeEfficiencyMean(float)));
    QObject::connect(model, SIGNAL(updateNumSearchers(int)), window, SLOT(setNumSearchers(int)));
    QObject::connect(model, SIGNAL(updateNumClusters(int)), window, SLOT(setNumClusters(int)));
    QObject::connect(model, SIGNAL(updateNumTargetsPerCluster(int)), window, SLOT(setNumTargetsPerCluster(int)));
    QObject::connect(model, SIGNAL(updateDistEfficiencyMean(float)), window, SLOT(setDistEfficiencyMean(float)));
    QObject::connect(model, SIGNAL(updateDistEfficiencyStd(float)), window, SLOT(setDistEfficiencyStd(float)));
    QObject::connect(model, SIGNAL(updateTimeEfficiencyStd(float)), window, SLOT(setTimeEfficiencyStd(float)));

    QObject::connect(model, SIGNAL(updateTimeResolution(float)), window, SLOT(setTimeResolution(float)));
    QObject::connect(model, SIGNAL(updateVolume(float)), window, SLOT(setVolume(float)));
    QObject::connect(model, SIGNAL(updateDensity(float)), window, SLOT(setDensity(float)));
    QObject::connect(model, SIGNAL(updateNumTargets(int)), window, SLOT(setNumTargets(int)));


    QObject::connect(window->total_contacts_value, SIGNAL(toggled(bool)), model, SLOT(setTargetsWithReplacement(bool)));
    window->total_contacts_value->toggle();


    QObject::connect(model, SIGNAL(updateClusterRadius(float)), window, SLOT(setClusterRadius(float)));
    QObject::connect(model, SIGNAL(updateProgressBar(int)), window, SLOT(setProgressBarPercentage(int)));
    QObject::connect(model, SIGNAL(updateHistogramValues(MyArray)), window, SLOT(setHistogramValues(MyArray)));


    QObject::connect(model, SIGNAL(updateTimeEffSample(FloatVector)), window, SLOT(setTimeEfficiencySample(FloatVector)));
    QObject::connect(model, SIGNAL(updateDistEffSample(FloatVector)), window, SLOT(setDistanceEfficiencySample(FloatVector)));
    QObject::connect(model, SIGNAL(updateHistogramMaxValue(float)), window, SLOT(setHistogramMaxValue(float)));

    QObject::connect(model, SIGNAL(addInputFilesToProcess(QStringList)), window, SLOT(addInputFilesToProcess(QStringList)));

    QObject::connect(model, SIGNAL(updateInputFileProcessed(int)), window, SLOT(setInputFileProcessed(int)));
    QObject::connect(model, SIGNAL(updateInputFileProcessing(int)), window, SLOT(setInputFileProcessing(int)));
    QObject::connect(model, SIGNAL(clearInputFilesToProcess()), window, SLOT(clearInputFilesToProcess()));

    QObject::connect(model, SIGNAL(updateDistExpended(float)), window, SLOT(setDistExpended(float)));
    QObject::connect(model, SIGNAL(updateTimeExpended(float)), window, SLOT(setTimeExpended(float)));

    //QObject::connect(visualizerThread, &QThread::started, vis, &Visualizer::updateGL);

    QObject::connect(vis, SIGNAL(StopModel()), model, SLOT(stop()));
    QObject::connect(vis, SIGNAL(StartModel()), model, SLOT(start()));

    //QObject::connect(vis, SIGNAL(generate()), model, SLOT(generateSearch()));
    QObject::connect(vis, SIGNAL(generate()), model, SLOT(analyse()));
    QObject::connect(vis, SIGNAL(VisPlaceSearchers()), model, SLOT(placeSearchers()));
    QObject::connect(vis, SIGNAL(VisSetSearchType(int)), model, SLOT(setSearchType(int)));
    QObject::connect(vis, SIGNAL(VisSetClusterRadius(float)), model, SLOT(setClusterRadius(float)));
    QObject::connect(vis, SIGNAL(VisFinishedPaintGL()), model, SLOT(paintingGLFinished()));


    window->setDetectionRadius(model->getDetectionRadius());
    window->setNumClusters(model->getNumClusters());
    window->setNumTargetsPerCluster(model->getNumTargetsPerCluster());
    window->setDisplayRadius(vis->getDisplayRadius());
    window->setClusterRadius(model->getClusterRadius());
   // window.setDisplayRadius(model->getNumSearchers());

    // Center the window
    //visualizerThread->start();

    string title  = "Power Search 3D: Version 8.0 Alpha";
    window->setWindowTitle(tr(title.c_str()));

    window->show();

    int width = qApp->desktop()->availableGeometry().width();
    int height = qApp->desktop()->availableGeometry().height();

    window->setGeometry(width*0.1, height*0.1, width-width*0.2, height-height*0.2);



}

void PowerSearch::quit()
{
    cout << "Application powersearch quitting" << endl;
    delete model;
    delete vis;
    delete window;
}
