#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QtWidgets>
#include <QWidget>
#include "Histogram.h"
#include "QCustomPlot.h"
#include <QTime>

using namespace std;

typedef QVector<int> MyArray;
typedef QVector<float> FloatVector; // typedef because of QT bug


QT_BEGIN_NAMESPACE
class QSlider;
class QLabel;
QT_END_NAMESPACE

class Visualizer;
class Model;

class MainWindow : public QWidget
{
    Q_OBJECT

public:
    MainWindow();
    ~MainWindow();
    Visualizer* getVisualizer();
    void setVisualizer(Visualizer* v);
 //   void setTimeEfficiencyValue( float v );
    void setDisplayRadius( float v );
    void setDetectionRadius( float v );
   // QCheckBox* total_contacts_value;
    QCheckBox* outliers_value, *dist_log_yaxis_value, *time_log_yaxis_value;
   // QCheckBox* unique_contacts_value;
    QRadioButton *destructive_radio, *nondestructive_radio, *unique_radio;
    QRadioButton *time_limit_radio, *dist_limit_radio;

    QTextEdit* myTextEdit;
    Model* getModel();
    QComboBox* search_type_combobox;

public slots:
    void setTimeExpended(float v);
    void setDistExpended(float v);
    void setTimeEfficiencyMean(float v);
    void setDistEfficiencyMean(float v);
    void setTimeEfficiencyStd(float v);
    void setDistEfficiencyStd(float v);
    void setDisplayOutliers(bool v);
    void setTimeResolution(float v);
    void setVolume(float v);
    void setDensity(float v);
    void setNumTargets(int v);
    void setNumSearchers(int v);
    void setNumClusters( int v );
    void setNumTargetsPerCluster( int v );
    void setDistScaleType(bool v);
    void setTimeScaleType(bool v);

    void setOutputDirPath(QString path);

    void setProgressBarPercentage(float);
    void setClusterRadius( float v );
    void setSpeedHistogramValues(MyArray values);
    void setSpeedHistogramMaxValue(float v);
    void setAngleHistogramValues(MyArray values);
    void setAngleHistogramMaxValue(float v);
    void addInputFilesToProcess(QStringList v);
    void setInputFileProcessed(int);
    void setInputFileProcessing(int);
    void clearInputFilesToProcess();
    void setTimeUnits(QString v);
    void setSpatialUnits(QString v);
    QStringList OpenFile();
    void addDataset();
    void modelFinished();

    void save();

    void updateTime();

    // Box plot slots
    void setTimeEfficiencySample(FloatVector values);
    void setDistanceEfficiencySample(FloatVector values);
    void setFirstContactSample(FloatVector values);

    signals:
    void setModelSavePath(QString);

protected:
    void keyPressEvent(QKeyEvent *event);

private:
    bool isNewDistSetID(int setid);
    bool isNewTimeSetID(int setid);
    bool isNewFirstContactSetID(int setid);


    vector<int> old_distsetids;
    vector<int> old_timesetids;
    vector<int> old_first_contact_setids;


    void setupBoxPlot();

    QSlider *createSlider();

    Visualizer* vis;
    QSlider* xSlider;
    QSlider* ySlider;
    QSlider* zSlider;

    QLineEdit* targets_detect_radius_value;
    QLineEdit* targets_display_radius_value;
    QLabel* time_efficiency_value;
    QLabel* dist_efficiency_value;
    QLabel* time_efficiency_std_value;
    QLabel* dist_efficiency_std_value;

    QLabel* time_expended_value;
    QLabel* dist_expended_value;
    QLineEdit* n_targets_per_cluster_value;
    QLineEdit* n_clusters_value;
    QLineEdit* n_searchers_value;
    QLineEdit* cluster_radius_value;
    QProgressBar* bar;
    Histogram* speed_histogram;
    Histogram* angle_histogram;
    vector<QListWidget*> input_files_to_process_lists;
    QLineEdit* time_resolution_value;
    QLabel* volume_value;
    QLabel* n_targets_value;
    QLabel* dist_efficiency_units_label;
    QLabel* time_efficiency_units_label;
    QLabel* density_value;
    QLineEdit* time_units_value;
    QLineEdit* spatial_units_value;
    QLineEdit* output_dir_path_value;



    QCustomPlot* timeCustomPlot;
    QCustomPlot* distCustomPlot;
    QCustomPlot* firstContactCustomPlot;


    vector<QCPStatisticalBox*> time_eff_samples;
    vector<QCPStatisticalBox*> dist_eff_samples;
    vector<QCPStatisticalBox*> first_contact_samples;


    QGridLayout* toolLayout;
    QVBoxLayout* queued_files_layout;

    int current_file_index;
    QTime start_time;
    QTime elapsed_time;
    int est_remaining_time_secs;

    QLabel* elapsed_time_value;
    QLabel* est_remaining_time_value;
    QTimer* timer;

    vector<QString> dataset_labels;
    int n_datasets;


    QVector<QBrush> boxBrushes;
    QString default_save_path;

    float density;
    float time_expended;
    float dist_expended;
    float time_resolution;

    bool display_outliers, dist_log_yaxis, time_log_yaxis;


    QVector<QString> time_x_axis_tick_vector_label;
    QVector<double> time_x_axis_tick_vector;
    QVector<QString> dist_x_axis_tick_vector_label;
    QVector<double> dist_x_axis_tick_vector;
    QVector<QString> first_contact_x_axis_tick_vector_label;
    QVector<double> first_contact_x_axis_tick_vector;


};
//! [0]

#endif
