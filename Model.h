#ifndef MODEL_H
#define MODEL_H

#include <QTimer>
#include "MainWindow.h"
#include <string>
#include <QVector>
#include <QMutex>
////#inlcude

#define PI 3.14159265359

using namespace std;

class Coordinate;
class Worker;
class Sphere;
class RectangularPrism;
class SearchSpace;
class MainWindow;

typedef QVector<int> MyArray;

class Model : public QObject
{

    Q_OBJECT

public:
  Model(float spr, float cr, int nc, int ntc, int sp, bool t, string app_path);
  SearchSpace* getSearchSpace();

  void writeIndividualTimeEfficiencies();
  void writeIndividualDistEfficiencies();
  void ParseTracks(string filename);
  void addSearchers(QString file_path);
  void clearPaths();
  ~Model();
 int calcTargetsFound();
 float getDetectionRadius();
 int getNumClusters();
 float getClusterRadius();
 int getNumTargetsPerCluster();
 int getNumSearchers();
 vector<Coordinate*> getHopkinsPoints();
 int getNumTargets();
 void removeTargets();
 void removeSearchers();
 void placeTargets();
 void generate();
 void run();
 int getSearchType();
 void doAnalysis();
 void saveSettings();
 void readSettings();
 void makeHistograms();
 // Output for ANOVA analysis
 void writeDistEfficienciesToFile(string file_path);
 void writeTimeEfficienciesToFile(string file_path);
 vector<double> calcOutliers(vector<double> values, double inner_fence, double median);
 string getWorkingPath();
 string getSavePath();
     bool gui_ready;
     // Two sample location Student's T-test
     double studentsTtest(vector<double> sample1, vector<double> sample2);
     double pValue(vector<double> );
     double mean( vector<double> values );
     double populationStandardDeviation( vector<double> values );
     double unbiasedSampleStandardDeviation( vector<double> values );
     double unbiasedSampleVariance( vector<double> values );
     double populationVariance( vector<double> values );
     double sum( vector<double> values );
     double median(vector<double> values);
     double lowerQuartile(vector<double> values);
     double upperQuartile(vector<double> values);
     double min(vector<double> values);
     double max(vector<double> values);
     double hopkins(vector<Coordinate*> locations);
     double mannwhitney(vector<double> values1, vector<double> values2);
     void studentsTtestDifferingSTDUsingBoost( string S1_label, double Sm1, double Sd1, unsigned Sn1, string S2_label, double Sm2, double Sd2, unsigned Sn2, double alpha);
     void studentsTtestUsingBoost( string S1_label, double Sm1, double Sd1, unsigned Sn1, string S2_label, double Sm2, double Sd2, unsigned Sn2, double alpha);
     pair<double, double> unbiasedSampleConfidenceInterval( vector<double> sample );
     pair<double, double> populationConfidenceInterval( vector<double> population );
     void printVectorContents( vector<double> values );
     void writeVectorContentsToCSVFile( vector<double> values, string filepath );
     void queueInputFiles(QStringList l);
     void lock();
     void unlock();
     void writeIndividualTrackEfficienciesToCSVFile( string filepath );
     QVector<float> packageSampleForHistogram(vector<double> data, int setid, int id);

     void addInputSet(QStringList* input_file_list);
     MainWindow* window;


 bool drawing;

signals:
  void updateNeeded();
  void updateTimeEfficiencyMean(float v);
  void updateDistEfficiencyMean(float v);
  void updateTimeEfficiencyStd(float v);
  void updateDistEfficiencyStd(float v);


  void updateTimeExpended(float v);
  void updateDistExpended(float v);
  void updateDetectRadius(float v);
  void updateNumSearchers(int v);
  void updateNumClusters(int v);
  void updateNumTargetsPerCluster(int v);
  void updateTimeResolution(float v);
  void updateClusterRadius(float v);
  void parseError(QString);
  void startWorker();
  void wakeupWorker();
  void updateProgressBar(float);
  void updateSpeedHistogramValues(MyArray values);
  void updateAngleHistogramValues(MyArray values);
  void updateDistEffSample(FloatVector values);
  void updateTimeEffSample(FloatVector values);
  void updateFirstContactSample(FloatVector values);
  void updateSpeedHistogramMaxValue(float v);
  void updateAngleHistogramMaxValue(float v);
  void addInputFilesToProcess(QStringList l);
  void updateInputFileProcessed(int);
  void updateInputFileProcessing(int);
  void clearInputFilesToProcess();
  void updateVolume(float);
  void updateDensity(float);
  void updateNumTargets(int);
  void modelFinished();

  public slots:
  void stop();
  void start();
  void reset();
  void generateSearch();
  void analyse();
  void repeat();
  void handleWorkerFinish();
  void placeSearchers();
  void setSearchType(int);
   void setClusterRadius(FloatVector v);
   void paintingGLFinished();
   void setTargetsWithReplacement(bool);
   void setUniqueTargets(bool);
   void setTimeLimited(bool);
   void setDistLimited(bool);
   void setSavePath(QString);

private:
  int max_targets;
  RectangularPrism* space;
  float max_time; // Maximum time searches are allowed (summed over all searchers)
  float max_distance; // Maximum distance searchers are allowed
  int time;
  float detection_radius;
  FloatVector cluster_radii;
  float cluster_radius;
  int n_clusters;
  int targets_per_cluster;
  int n_targets_found;
  float distance_used;
  float time_used;
  int search_pattern;
  float time_resolution;

  bool targets_with_replacement;
  bool unique_targets;

//  bool time_efficiency;
    QTimer* iTimer;
    bool stopped;
    QThread* workerThread;
    Worker* worker;
    int num_targets_found;
    bool safe_to_update_positions;
    float dist_efficiency;
    float time_efficiency;
    float first_contact;
    string save_path;
    string working_path;
    string app_path;
    QVector<int> velocity_distribution;

    vector<QStringList*> input_sets;

    vector<Coordinate*> hopkins_points;

    float reference_volume;
    float scale_number_of_targets;
    bool idealized_search;

    vector< vector<double>* > mean_speeds;
    vector< vector<double>* > distances_used;
    vector< vector<double>* > times_used;
    vector< vector<double>* > targets_found;
    vector< vector<double>* > dist_efficiencies;
    vector< vector<double>* > first_contacts;
    vector< vector<double>* > time_efficiencies;
    vector<double> time_efficiencies_replica_means;
    vector<double> dist_efficiencies_replica_means;
    vector<double> first_contact_replica_means;
    vector< vector<double>* > hopkins_values;
    QStringList file_factors; // Statistical factors read from input file
    QStringList cluster_factors; // Statistical factors including cluster size

    QStringList dataset_labels;

    bool limit_by_time;

    QMutex mutex;

};

#endif // MODEL_H
