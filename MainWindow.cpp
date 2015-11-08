
#include "Model.h"
#include "Visualizer.h"
#include "MainWindow.h"
#include <string>
#include <boost/lexical_cast.hpp>
#include "Histogram.h"
#include "QDebugStream.h"

using namespace std;
using namespace boost;

//! [0]
MainWindow::MainWindow()
{
    #if defined(__APPLE__)
    default_save_path += "/MyFile";
    #endif

    est_remaining_time_secs = 0;
    display_outliers = false;
    dist_log_yaxis = false;
    dist_log_yaxis = false;

    density = 0;
    time_expended = 0;
    dist_expended = 0;
    time_resolution = 0;

    default_save_path = "";

    n_datasets = 0;

    setWindowIcon(QIcon("powersearch.icns"));

    vis = new Visualizer();

    timer = new QTimer();

    connect(timer, SIGNAL(timeout()), this, SLOT(updateTime()));
    timer->start(1000);

    xSlider = createSlider();
    ySlider = createSlider();
    zSlider = createSlider();

    connect(xSlider, SIGNAL(valueChanged(int)), vis, SLOT(setXRotation(int)));
    connect(vis, SIGNAL(xRotationChanged(int)), xSlider, SLOT(setValue(int)));
    connect(ySlider, SIGNAL(valueChanged(int)), vis, SLOT(setYRotation(int)));
    connect(vis, SIGNAL(yRotationChanged(int)), ySlider, SLOT(setValue(int)));
    connect(zSlider, SIGNAL(valueChanged(int)), vis, SLOT(setZRotation(int)));
    connect(vis, SIGNAL(zRotationChanged(int)), zSlider, SLOT(setValue(int)));
//! [0]

//! [1]


    QGroupBox *searchTypeGroupBox = new QGroupBox(tr("Contact Properties"));

         destructive_radio = new QRadioButton(tr("Destructive"));
         nondestructive_radio = new QRadioButton(tr("Non-destructive"));
         unique_radio = new QRadioButton(tr("Unique"));

         QVBoxLayout *searchTypevbox = new QVBoxLayout;
         searchTypevbox->addWidget(destructive_radio);
         searchTypevbox->addWidget(nondestructive_radio);
         searchTypevbox->addWidget(unique_radio);
         searchTypevbox->addStretch(1);
         searchTypeGroupBox->setLayout(searchTypevbox);


         QGroupBox *limitTypeGroupBox = new QGroupBox(tr("Limit Properties"));

              time_limit_radio = new QRadioButton(tr("Time"));
              dist_limit_radio = new QRadioButton(tr("Distance"));

              QVBoxLayout *limitTypevbox = new QVBoxLayout;
              limitTypevbox->addWidget(time_limit_radio);
              limitTypevbox->addWidget(dist_limit_radio);
              limitTypevbox->addStretch(1);
              limitTypeGroupBox->setLayout(limitTypevbox);


    myTextEdit = new QTextEdit(this);

    QVBoxLayout* mainLayout = new QVBoxLayout;
    QHBoxLayout* visualization_and_tools_layout = new QHBoxLayout;
    QHBoxLayout* plots_layout = new QHBoxLayout;
    QHBoxLayout* visualizations_layout = new QHBoxLayout;


    speed_histogram = new Histogram(this);
    angle_histogram = new Histogram(this);


    //histogram->setBins(bins);

    speed_histogram->show();
    angle_histogram->show();

    toolLayout = new QGridLayout;

    /*
    QLabel* total_contacts_label = new QLabel("Contacts With Replacement");
    total_contacts_value = new QCheckBox();

    QLabel* unique_contacts_label = new QLabel("Unique Contacts");
    unique_contacts_value = new QCheckBox();
    */

    QLabel* outliers_label = new QLabel("Show Outliers");
    outliers_value = new QCheckBox();
    QObject::connect(outliers_value, SIGNAL(toggled(bool)), this, SLOT(setDisplayOutliers(bool)));

    QLabel* dist_log_yaxis_label = new QLabel("Log Distance Y-Axis");
    dist_log_yaxis_value = new QCheckBox();
    QObject::connect(dist_log_yaxis_value, SIGNAL(toggled(bool)), this, SLOT(setDistScaleType(bool)));

    QLabel* time_log_yaxis_label = new QLabel("Log Time Y-Axis");
    time_log_yaxis_value = new QCheckBox();
    QObject::connect(time_log_yaxis_value, SIGNAL(toggled(bool)), this, SLOT(setTimeScaleType(bool)));


    QLabel* volume_label = new QLabel("Search Space Volume");
    volume_value = new QLabel("0");

    QLabel* density_label = new QLabel("Target Density");
    density_value = new QLabel("0");

    QLabel* n_targets_label = new QLabel("Number of Targets");
    n_targets_value = new QLabel("0");

    QLabel* time_resolution_label = new QLabel("Time Resolution");
    time_resolution_value = new QLineEdit("0");

    QLabel* n_searchers_label = new QLabel("Number of Searchers");
    n_searchers_value = new QLineEdit("0");

    QLabel* n_clusters_label = new QLabel("Number of Target Clusters");
    n_clusters_value = new QLineEdit("Value");


    QLabel* n_targets_per_cluster_label = new QLabel("Number of Targets per Cluster");
    n_targets_per_cluster_value = new QLineEdit("Value");


    QLabel* cluster_radius_label = new QLabel("Cluster Radius");
    cluster_radius_value = new QLineEdit("Value");

    connect(cluster_radius_value, SIGNAL(textChanged(const QString &)), vis, SLOT(GUISetClusterRadius(const QString &)));
    QLabel* targets_detect_radius_label = new QLabel("Target Detect Radius");
    targets_detect_radius_value = new QLineEdit("0");

    connect(targets_detect_radius_value, SIGNAL(textChanged(const QString &)), vis, SLOT(GUISetTargetDetectRadius(const QString &)));

    QLabel* targets_display_radius_label = new QLabel("Target Display Radius");

    targets_display_radius_value = new QLineEdit("0");

    connect(targets_display_radius_value, SIGNAL(textChanged(const QString &)), vis, SLOT(GUISetTargetDisplayRadius(const QString &)));

    QLabel* spatial_units_label = new QLabel("Spatial Units");
    spatial_units_value = new QLineEdit("\u00b5m");
    connect(spatial_units_value, SIGNAL(textChanged(const QString &)), this, SLOT(setSpatialUnits(const QString &)));

    QLabel* time_units_label = new QLabel("Time Units");
    time_units_value = new QLineEdit("s");
    connect(time_units_value, SIGNAL(textChanged(const QString &)), this, SLOT(setTimeUnits(const QString &)));

    QLabel* time_expended_label = new QLabel("Time Expended");

    time_expended_value = new QLabel("0");

    QLabel* time_expended_units_label = new QLabel("s");

    QLabel* dist_expended_label = new QLabel("Distance Expended");
    dist_expended_value = new QLabel("0");

    QLabel* dist_expended_units_label = new QLabel(spatial_units_value->text());

    time_efficiency_std_value = new QLabel("0");
    dist_efficiency_std_value = new QLabel("0");


    QLabel* time_efficiency_label = new QLabel("Time Efficiency");
    time_efficiency_value = new QLabel("0");
    time_efficiency_units_label = new QLabel(QString::fromUtf8("targets/")+time_units_value->text());


    QLabel* dist_efficiency_label = new QLabel("Distance Efficiency");
    dist_efficiency_value = new QLabel("0");
    dist_efficiency_units_label = new QLabel(QString::fromUtf8("targets/")+spatial_units_value->text());

    QPushButton* place_searchers_button = new QPushButton("Place Searchers");

    connect(place_searchers_button, SIGNAL(clicked()), vis, SLOT(GUIPlaceSearchers()));

    QLabel* output_dir_path_label = new QLabel("Path to Output Directory");
    output_dir_path_value = new QLineEdit(default_save_path);
    connect(output_dir_path_value, SIGNAL(textChanged(const QString &)), this, SLOT(setOutputDirPath(const QString &)));


    int tools_row = 0;

    toolLayout->addWidget(n_clusters_label, tools_row, 1);
    toolLayout->addWidget(n_clusters_value, tools_row, 2);

    toolLayout->addWidget(searchTypeGroupBox, tools_row, 3, 3, 3);

    //toolLayout->addWidget(total_contacts_label, tools_row, 3);
    //toolLayout->addWidget(total_contacts_value, tools_row, 4);

    tools_row++;

    toolLayout->addWidget(n_targets_per_cluster_label, tools_row, 1);
    toolLayout->addWidget(n_targets_per_cluster_value, tools_row, 2);

    //toolLayout->addWidget(unique_contacts_label, tools_row, 3);
    //toolLayout->addWidget(unique_contacts_value, tools_row, 4);

    tools_row++;

    toolLayout->addWidget(cluster_radius_label, tools_row, 1);
    toolLayout->addWidget(cluster_radius_value, tools_row, 2);


    tools_row++;


    toolLayout->addWidget(targets_display_radius_label, tools_row, 1);
    toolLayout->addWidget(targets_display_radius_value, tools_row, 2);

    toolLayout->addWidget(outliers_label, tools_row, 3);
    toolLayout->addWidget(outliers_value, tools_row, 4);

    tools_row++;

    toolLayout->addWidget(spatial_units_label, tools_row, 1);

    toolLayout->addWidget(spatial_units_value, tools_row, 2);


    toolLayout->addWidget(time_log_yaxis_label, tools_row, 3);

    toolLayout->addWidget(time_log_yaxis_value, tools_row, 4);


    tools_row++;


    toolLayout->addWidget(time_units_label, tools_row, 1);
    toolLayout->addWidget(time_units_value, tools_row, 2);

    toolLayout->addWidget(dist_log_yaxis_label, tools_row, 3);
    toolLayout->addWidget(dist_log_yaxis_value, tools_row, 4);


    tools_row++;


    toolLayout->addWidget(targets_detect_radius_label, tools_row, 1);
    toolLayout->addWidget(targets_detect_radius_value, tools_row, 2);

    toolLayout->addWidget(limitTypeGroupBox, tools_row, 3, 3, 3);

    tools_row++;

    toolLayout->addWidget(n_searchers_label, tools_row, 1);
    toolLayout->addWidget(n_searchers_value, tools_row, 2);

    //toolLayout->addWidget(place_searchers_button, tools_row, 3);



    tools_row++;

    toolLayout->addWidget(time_resolution_label, tools_row, 1);
    toolLayout->addWidget(time_resolution_value, tools_row, 2);

    tools_row++;

    toolLayout->addWidget(volume_label, tools_row, 1);
    toolLayout->addWidget(volume_value, tools_row, 2);

    tools_row++;

    toolLayout->addWidget(density_label, tools_row, 1);
    toolLayout->addWidget(density_value, tools_row, 2);

    tools_row++;

    toolLayout->addWidget(n_targets_label, tools_row, 1);
    toolLayout->addWidget(n_targets_value, tools_row, 2);

    tools_row++;

    toolLayout->addWidget(xSlider, tools_row, 1);
    toolLayout->addWidget(ySlider, tools_row, 2);
    toolLayout->addWidget(zSlider, tools_row, 3);

    xSlider->setValue(15 * 16);
    ySlider->setValue(345 * 16);
    zSlider->setValue(0 * 16);

    QPushButton* zoomin_button = new QPushButton("Zoom In");
    QPushButton* zoomout_button = new QPushButton("Zoom Out");
    QPushButton* openfile_button = new QPushButton("Add Dataset");

    QPushButton* panup_button = new QPushButton("Pan Up");
    QPushButton* pandown_button = new QPushButton("Pan Down");
    QPushButton* panleft_button = new QPushButton("Pan Left");
    QPushButton* panright_button = new QPushButton("Pan Right");


    connect(zoomin_button, SIGNAL(clicked()), vis, SLOT(ZoomIn()));
    connect(zoomout_button, SIGNAL(clicked()), vis, SLOT(ZoomOut()));
    connect(openfile_button, SIGNAL(clicked()), this, SLOT(addDataset()));

    connect(panup_button, SIGNAL(clicked()), vis, SLOT(PanUp()));
    connect(pandown_button, SIGNAL(clicked()), vis, SLOT(PanDown()));
    connect(panleft_button, SIGNAL(clicked()), vis, SLOT(PanLeft()));
    connect(panright_button, SIGNAL(clicked()), vis, SLOT(PanRight()));

    tools_row++;

    toolLayout->addWidget(zoomin_button, tools_row, 1);
    toolLayout->addWidget(zoomout_button, tools_row, 2);
    toolLayout->addWidget(openfile_button, tools_row, 3);

    tools_row++;

    toolLayout->addWidget(panup_button, tools_row, 2);

    tools_row++;

    toolLayout->addWidget(panleft_button, tools_row, 1);

    toolLayout->addWidget(panright_button, tools_row, 3);

    tools_row++;

    toolLayout->addWidget(pandown_button, tools_row, 2);

    QPushButton* generate_button = new QPushButton("Analyse");
    QPushButton* save_button = new QPushButton("Set Output Path");
    search_type_combobox = new QComboBox();
    search_type_combobox->addItem("Observed Motion");
    search_type_combobox->addItem("Brownian Motion");
    search_type_combobox->addItem("LogNormal PDF");
    search_type_combobox->addItem("Gen. Pareto PDF");
    search_type_combobox->addItem("Power Law PDF");
    search_type_combobox->addItem("Exponential PDF");
    search_type_combobox->addItem("Gamma PDF");
    search_type_combobox->addItem("CRW");
    search_type_combobox->addItem("LogMCRW");
    search_type_combobox->addItem("GammaMCRW");
    search_type_combobox->addItem("Sample Observed Speeds (Bootstrap)");
    search_type_combobox->addItem("Sample Observed Angles and Speeds (Bootstrap)");

    connect(save_button, SIGNAL(clicked()), this, SLOT(save()));
    connect(generate_button, SIGNAL(clicked()), vis, SLOT(GUIgenerate()));

    tools_row++;

    toolLayout->addWidget(generate_button, tools_row, 1);
    toolLayout->addWidget(search_type_combobox, tools_row, 2);
    toolLayout->addWidget(save_button, tools_row, 3);

    tools_row++;

    toolLayout->addWidget(time_expended_label, tools_row, 1);
    toolLayout->addWidget(time_expended_value, tools_row, 2);
    toolLayout->addWidget(output_dir_path_label, tools_row, 3);

    // toolLayout->addWidget(time_expended_units_label, tools_row, 3);

    tools_row++;

    toolLayout->addWidget(dist_expended_label, tools_row, 1);
    toolLayout->addWidget(dist_expended_value, tools_row, 2);
    toolLayout->addWidget(output_dir_path_value, tools_row, 3, 1, 2);

 //   toolLayout->addWidget(dist_expended_units_label, tools_row, 3);
    setDistExpended(0);
    setTimeExpended(0);

    tools_row++;

    toolLayout->addWidget(time_efficiency_label, tools_row, 1);
    toolLayout->addWidget(time_efficiency_value, tools_row, 2);
    toolLayout->addWidget(time_efficiency_std_value, tools_row, 3);
    toolLayout->addWidget(time_efficiency_units_label, tools_row, 4);

    tools_row++;

    toolLayout->addWidget(dist_efficiency_label, tools_row, 1);
    toolLayout->addWidget(dist_efficiency_value, tools_row, 2);
    toolLayout->addWidget(dist_efficiency_std_value, tools_row, 3);
    toolLayout->addWidget(dist_efficiency_units_label, tools_row, 4);

    bar = new QProgressBar(0);
    bar->setRange(0, 100);
    bar->setValue(0);
    bar->show();

    speed_histogram->toggled();
    speed_histogram->setMinimumHeight(100);
    speed_histogram->setXUnits(spatial_units_value->text() + "/" + time_units_value->text());

    angle_histogram->toggled();
    angle_histogram->setMinimumHeight(100);
    angle_histogram->setXUnits("\u03C0");

    tools_row++;

    toolLayout->addWidget(new QLabel("Speed Histogram"), tools_row, 1, 1, 1);
    toolLayout->addWidget(new QLabel("Angle Histogram"), tools_row, 3, 1, 1);

    tools_row++;

    toolLayout->addWidget(speed_histogram, tools_row, 1, 6, 2);
    toolLayout->addWidget(angle_histogram, tools_row, 3, 6, 2);

    tools_row+=6;

    toolLayout->addWidget(bar, tools_row, 1);


    //
    // Begin Time boxplot
    //

    timeCustomPlot= new QCustomPlot();
    distCustomPlot= new QCustomPlot();
    firstContactCustomPlot= new QCustomPlot();

    timeCustomPlot->axisRect()->setRangeZoomFactor(0.99);
    distCustomPlot->axisRect()->setRangeZoomFactor(0.99);
    firstContactCustomPlot->axisRect()->setRangeZoomFactor(0.99);

    // Predefined brushes
    QBrush boxBrush1(QColor(60, 60, 255, 100));
    boxBrush1.setStyle(Qt::NoBrush); // make it look oldschool

    QBrush boxBrush2(QColor(60, 255, 60, 100));
    boxBrush2.setStyle(Qt::SolidPattern); // make it look oldschool

    QBrush boxBrush3(QColor(255, 60, 60, 100));
    boxBrush3.setStyle(Qt::FDiagPattern); // make it look oldschool

    QBrush boxBrush4(QColor(60, 100, 255, 100));
    boxBrush4.setStyle(Qt::CrossPattern); // make it look oldschool

    QBrush boxBrush5(QColor(100, 255, 60, 100));
    boxBrush5.setStyle(Qt::BDiagPattern); // make it look oldschool

    QBrush boxBrush6(QColor(255, 60, 150, 100));
    boxBrush6.setStyle(Qt::VerPattern); // make it look oldschool

    QBrush boxBrush7(QColor(255, 0, 200, 100));
    boxBrush7.setStyle(Qt::DiagCrossPattern); // make it look oldschool

    QBrush boxBrush8(QColor(0, 200, 60, 100));
    boxBrush8.setStyle(Qt::HorPattern); // make it look oldschool

    boxBrushes.push_back(boxBrush1);
    boxBrushes.push_back(boxBrush2);
    boxBrushes.push_back(boxBrush3);
    boxBrushes.push_back(boxBrush4);
    boxBrushes.push_back(boxBrush5);
    boxBrushes.push_back(boxBrush6);
    boxBrushes.push_back(boxBrush7);
    boxBrushes.push_back(boxBrush8);

    // create empty statistical box plottables:

    /*

    // One extra at each side for padding on axis rescale
    for (int i = 0; i < 17; i++)
         time_eff_samples.push_back(new QCPStatisticalBox(timeCustomPlot->xAxis, timeCustomPlot->yAxis));

    QBrush observed_boxBrush(QColor(60, 60, 255, 100));
    observed_boxBrush.setStyle(Qt::DiagCrossPattern); // make it look oldschool

    QBrush brownian_boxBrush(QColor(60, 255, 60, 100));
    brownian_boxBrush.setStyle(Qt::Dense6Pattern); // make it look oldschool

    QBrush lognormal_boxBrush(QColor(255, 60, 60, 100));
    lognormal_boxBrush.setStyle(Qt::FDiagPattern); // make it look oldschool

    timeCustomPlot->legend->setVisible(true);

    QFont legendFont = font();  // start out with MainWindow's font..
    //legendFont.setPointSize(9); // and make a bit smaller for legend
    timeCustomPlot->legend->setFont(legendFont);


       for (int i = 1; i < 6; i++)
       {
           time_eff_samples[i]->setOutlierStyle(QCPScatterStyle(QCPScatterStyle::ssPlus, QColor(60, 60, 255, 100), 5));
           time_eff_samples[i]->setMeanStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::black, 8));

           time_eff_samples[i]->setBrush(observed_boxBrush);
        }

       for (int i = 6; i < 11; i++)
       {
           time_eff_samples[i]->setOutlierStyle(QCPScatterStyle(QCPScatterStyle::ssPlus, QColor(60, 255, 60, 100), 5));
           time_eff_samples[i]->setMeanStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::black, 8));
           time_eff_samples[i]->setBrush(brownian_boxBrush);
        }

       for (int i = 11; i < 16; i++)
       {
           time_eff_samples[i]->setOutlierStyle(QCPScatterStyle(QCPScatterStyle::ssPlus, QColor(255, 60, 60, 100), 5));
           time_eff_samples[i]->setMeanStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::black, 8));
           time_eff_samples[i]->setBrush(lognormal_boxBrush);
        }


       QVector<double> x_axis_tick_vector;

    for (int i = 0; i < time_eff_samples.size(); i++)
    {
        timeCustomPlot->addPlottable(time_eff_samples[i]);
        x_axis_tick_vector << i;

// Have to initialize the time_eff_samples here (not sure why)
        time_eff_samples[i]->setKey(i);
        time_eff_samples[i]->setMedian(0.0);
        time_eff_samples[i]->setLowerQuartile(0.0);
        time_eff_samples[i]->setUpperQuartile(0.0);
        time_eff_samples[i]->setMinimum(0.0);
        time_eff_samples[i]->setMaximum(0.0);
    }

    x_axis_tick_vector_label << "" << "10" << "20" << "30" << "40" << "50" << "10" << "20" << "30" << "40" << "50" << "10" << "20" << "30" << "40" << "50";


    // prepare manual x axis labels:
    timeCustomPlot->xAxis->setSubTickCount(0);
    timeCustomPlot->xAxis->setTickLength(0, 4);
    timeCustomPlot->xAxis->setTickLabelRotation(0);
    timeCustomPlot->xAxis->setAutoTicks(false);
    timeCustomPlot->xAxis->setAutoTickLabels(false);
    timeCustomPlot->xAxis->setTickVector(x_axis_tick_vector);
    timeCustomPlot->xAxis->setTickVectorLabels(x_axis_tick_vector_label);

    // prepare axes:
    timeCustomPlot->yAxis->setLabel("Time Efficiency (Targets/"+time_units_value->text()+")");
    timeCustomPlot->rescaleAxes();
    timeCustomPlot->xAxis->scaleRange(1.7, timeCustomPlot->xAxis->range().center());
    timeCustomPlot->yAxis->setRange(0, 7.0);
    timeCustomPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    QString s = "Target Cluster Radius ("+spatial_units_value->text()+")\n(Hopkins Statistic)";
    timeCustomPlot->xAxis->setLabel(QString::fromUtf8(s.toStdString().c_str()));


*/

    //
    // End boxplot
    //

    //timeCustomPlot->setMinimumWidth(600);
   // timeCustomPlot->setMinimumHeight(600);

    //
    // Begin Distance Effciency boxplot
    //



    // create empty statistical box plottables:
/*
    for (int i = 0; i < 17; i++)
        dist_eff_samples.push_back(new QCPStatisticalBox(distCustomPlot->xAxis, distCustomPlot->yAxis));

       for (int i = 1; i < 6; i++)
       {
          dist_eff_samples[i]->setOutlierStyle(QCPScatterStyle(QCPScatterStyle::ssPlus, QColor(60, 60, 255, 100), 5));
          dist_eff_samples[i]->setMeanStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::black, 8));

            dist_eff_samples[i]->setBrush(observed_boxBrush);
        }

       for (int i = 6; i < 11; i++)
       {
           dist_eff_samples[i]->setOutlierStyle(QCPScatterStyle(QCPScatterStyle::ssPlus,QColor(60, 255, 60, 100), 5));
           dist_eff_samples[i]->setMeanStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::black, 8));
            dist_eff_samples[i]->setBrush(brownian_boxBrush);
}

       for (int i = 11; i < 16; i++)
       {
           dist_eff_samples[i]->setOutlierStyle(QCPScatterStyle(QCPScatterStyle::ssPlus, QColor(255, 60, 60, 100), 5));
           dist_eff_samples[i]->setMeanStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::black, 8));

            dist_eff_samples[i]->setBrush(lognormal_boxBrush);
}

    for (int i = 0; i < dist_eff_samples.size(); i++)
    {
        distCustomPlot->addPlottable(dist_eff_samples[i]);
        x_axis_tick_vector << i;


// Have to initialize the dist_eff_samples here (not sure why)
        dist_eff_samples[i]->setKey(i);
        dist_eff_samples[i]->setMedian(0.0);
        dist_eff_samples[i]->setLowerQuartile(0.0);
        dist_eff_samples[i]->setUpperQuartile(0.0);
        dist_eff_samples[i]->setMinimum(0.0);
        dist_eff_samples[i]->setMaximum(0.0);
    }



    // prepare manual x axis labels:
    distCustomPlot->xAxis->setSubTickCount(0);
    distCustomPlot->xAxis->setTickLength(0, 4);
    distCustomPlot->xAxis->setTickLabelRotation(0);
    distCustomPlot->xAxis->setAutoTicks(false);
    distCustomPlot->xAxis->setAutoTickLabels(false);
    distCustomPlot->xAxis->setTickVector(x_axis_tick_vector);
    distCustomPlot->xAxis->setTickVectorLabels(x_axis_tick_vector_label);

    // prepare axes:
    s = "Distance Efficiency (Targets/"+spatial_units_value->text()+")";
    distCustomPlot->yAxis->setLabel(s);
    distCustomPlot->rescaleAxes();
    distCustomPlot->xAxis->scaleRange(1.7, distCustomPlot->xAxis->range().center());
    distCustomPlot->yAxis->setRange(0, 7.0);
    distCustomPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);

    s = "Target Cluster Radius ("+spatial_units_value->text()+")\n(Hopkins Statistic)";
    distCustomPlot->xAxis->setLabel(s);

    //
    // End boxplot
    //

    dataset_labels.push_back("One");
    dataset_labels.push_back("Two");
    dataset_labels.push_back("Three");

    // setup the efficiency plot legend
    timeCustomPlot->legend->removeItem(0);

    time_eff_samples[1]->setName(dataset_labels[0]);

    for (int i = 0; i < 4; i++ )
    {
        timeCustomPlot->legend->removeItem(1);
    }

    time_eff_samples[6]->setName(dataset_labels[1]);

    for (int i = 0; i < 4; i++ )
    {
        timeCustomPlot->legend->removeItem(2);
    }

    time_eff_samples[11]->setName(dataset_labels[2]);


    for (int i = 0; i < 8; i++ )
    {
        timeCustomPlot->legend->removeItem(3);
    }

    timeCustomPlot->replot();



   // distCustomPlot->setMinimumWidth(600);
   // distPlot->setMinimumHeight(600);

    timeCustomPlot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignRight|Qt::AlignBottom);


    distCustomPlot->legend->setVisible(true);
    distCustomPlot->legend->setFont(legendFont);
        // setup the efficiency plot legend
        distCustomPlot->legend->removeItem(0);


        for (int i = 0; i < 4; i++ )
        {
            distCustomPlot->legend->removeItem(1);
        }

        dist_eff_samples[6]->setName(dataset_labels[1]);

        for (int i = 0; i < 4; i++ )
        {
            distCustomPlot->legend->removeItem(2);
        }

        dist_eff_samples[11]->setName(dataset_labels[2]);


        for (int i = 0; i < 8; i++ )
        {
            distCustomPlot->legend->removeItem(3);
        }

        distCustomPlot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignRight|Qt::AlignBottom);


        distCustomPlot->replot();
*/

    QVBoxLayout* secondaryLayout = new QVBoxLayout;

    distCustomPlot->setMinimumHeight(300);
    timeCustomPlot->setMinimumHeight(300);
    firstContactCustomPlot->setMinimumHeight(300);

    visualization_and_tools_layout->addWidget(vis);
    visualization_and_tools_layout->addLayout(toolLayout);

    plots_layout->addWidget(distCustomPlot);
    plots_layout->addWidget(timeCustomPlot);
    plots_layout->addWidget(firstContactCustomPlot);

    mainLayout->addLayout(visualization_and_tools_layout);
    mainLayout->addLayout(plots_layout);

    tools_row++;

    toolLayout->addWidget(new QLabel("Time Elapsed: "), tools_row, 1);
    toolLayout->addWidget(new QLabel("Est. Time Remaining: "), tools_row, 3);

    elapsed_time_value = new QLabel();
    est_remaining_time_value = new QLabel();

    toolLayout->addWidget(elapsed_time_value, tools_row, 2);
    toolLayout->addWidget(est_remaining_time_value, tools_row, 4);

    tools_row++;

    queued_files_layout = new QVBoxLayout;
    //toolLayout->addLayout(queued_files_layout, tools_row, 1, 6, 6);
    visualization_and_tools_layout->addLayout(queued_files_layout);

    myTextEdit->setMinimumWidth(500);
    visualization_and_tools_layout->addWidget(myTextEdit);

    setLayout(mainLayout);

    QList<QWidget*> list = this->findChildren<QWidget *>();
    foreach(QWidget *w, list)
    {
       w->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Expanding);
    }

    list = toolLayout->findChildren<QWidget *>();
    foreach(QWidget *w, list)
    {
       w->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Expanding);
    }



}
//! [1]

//! [2]
QSlider *MainWindow::createSlider()
{
    QSlider *slider = new QSlider(Qt::Vertical);
    slider->setRange(0, 360 * 16);
    slider->setSingleStep(16);
    slider->setPageStep(15 * 16);
    slider->setTickInterval(15 * 16);
    slider->setTickPosition(QSlider::TicksRight);
    return slider;
}
//! [2]

void MainWindow::keyPressEvent(QKeyEvent *e)
{
    if (e->key() == Qt::Key_Escape)
        close();
    else
        QWidget::keyPressEvent(e);
}

Visualizer* MainWindow::getVisualizer()
{
    return vis;
}

void MainWindow::setVisualizer( Visualizer* v )
{
    vis = v;
}

void MainWindow::setTimeExpended( float v )
{
    //cout << "MainWindow::setTimeExpended() called" << endl;
    time_expended = v;
    QString text = QString::number(v, 'g', 3)+time_units_value->text();
    time_expended_value->setText(text);
}

void MainWindow::setDistExpended( float v )
{
    dist_expended  = v;
    cout << "MainWindow::setDistExpended() called with value " << v << endl;
    QString text = QString::number(v, 'g', 3)+spatial_units_value->text();
    dist_expended_value->setText(text);
}

void MainWindow::setTimeEfficiencyStd( float v )
{
    QString text = QString::number(v, 'g', 3);
    time_efficiency_std_value->setText(text);
}

void MainWindow::setDistEfficiencyStd( float v )
{
    QString text = QString::number(v, 'g', 3);
    dist_efficiency_std_value->setText(text);
}

void MainWindow::setTimeEfficiencyMean( float v )
{
    QString text = QString::number(v, 'g', 3);
    time_efficiency_value->setText(text);
}

void MainWindow::setClusterRadius( float v )
{
    QString text = QString::number(v, 'g', 3);
    cluster_radius_value->setText(text);
}

void MainWindow::setDistEfficiencyMean( float v )
{
    QString text = QString::number(v, 'g', 3);
    dist_efficiency_value->setText(text);
}

void MainWindow::setDisplayRadius(float v)
{
    QString targets_display_radius_str = QString::number(v, 'g', 5);
    targets_display_radius_value->setText( targets_display_radius_str );
}

void MainWindow::setDetectionRadius(float v)
{
    QString targets_detection_radius_str = QString::number(v, 'g', 5);
    targets_detect_radius_value->setText( targets_detection_radius_str );
}

void MainWindow::setNumClusters(int v)
{
    QString str = QString::number(v, 'g', 3);
    n_clusters_value->setText( str );
}

void MainWindow::setNumTargetsPerCluster(int v)
{
    QString str = QString::number(v, 'g', 3);
    n_targets_per_cluster_value->setText( str );
}

void MainWindow::setNumSearchers(int v)
{
    QString str = QString::number(v, 'g', 3);
    n_searchers_value->setText( str );
}

void MainWindow::setTimeResolution(float v)
{
    time_resolution = v;
    QString str = QString::number(v, 'g', 5);
    str += time_units_value->text();
    time_resolution_value->setText( str );
}

void MainWindow::setVolume(float v)
{
    QString str = QString::number(v, 'g', 3);
    str += " (Ref. Vol.: 6.3e+06) "+spatial_units_value->text()+"\u00b3";
    volume_value->setText( str );
}

void MainWindow::setDensity(float v)
{
    density = v;
    QString str = QString::number(v, 'g', 3);
    str += " targets/"+spatial_units_value->text()+"\u00b3";
    density_value->setText( str );
}

void MainWindow::setNumTargets(int v)
{
    QString str = QString::number(v, 'g', 3);
    n_targets_value->setText( str );
}

MainWindow::~MainWindow()
{
    cout << "MainWindow: cleaning up" << endl;
    cout << "MainWindow: Saving plots" << endl;

    stringstream time_ss, dist_ss, diff_ss;
    time_ss << default_save_path.toStdString() << "/time.pdf";
    dist_ss << default_save_path.toStdString() << "/dist.pdf";
    diff_ss << default_save_path.toStdString() << "/diff.pdf";

    QString time_path = QString::fromStdString(time_ss.str());
    QString dist_path = QString::fromStdString(dist_ss.str());
    QString diff_path = QString::fromStdString(diff_ss.str());

    timeCustomPlot->savePdf(time_path);
    distCustomPlot->savePdf(dist_path);
    firstContactCustomPlot->savePdf(diff_path);

    cout << "MainWindow: saved time efficiency plot to " << time_path.toStdString() << endl;
    cout << "MainWindow: saved distance efficiency plot to " << dist_path.toStdString() << endl;
    cout << "MainWindow: saved difference efficiency plot to " << diff_path.toStdString() << endl;
}

void MainWindow::setProgressBarPercentage(float percentage)
{
    //cout << "MainWindow: Progress Bar Updated: " << percentage << endl;
    bar->setValue(ceil(percentage));

     if (elapsed_time.elapsed() == 0) elapsed_time.start();
     double rate_in_seconds = ((elapsed_time.elapsed()/1000.0f))/percentage;

     //cout << "MainWindow::setProgressBarPercentage(): rate in seconds: " << rate_in_seconds << endl;

     est_remaining_time_secs = ((100-percentage)*rate_in_seconds);

     if (ceil(percentage) == 100) timer->stop();
}

void MainWindow::setSpeedHistogramMaxValue(float v)
{
    speed_histogram->setMaxValue(v);
}

void MainWindow::setSpeedHistogramValues(MyArray values)
{
    speed_histogram->setBins(values);

}

void MainWindow::setAngleHistogramMaxValue(float v)
{
    angle_histogram->setMaxValue(v);
}

void MainWindow::setAngleHistogramValues(MyArray values)
{
    angle_histogram->setBins(values);

}

void MainWindow::addInputFilesToProcess(QStringList list)
{

    QListWidget* input_files_to_process_list = new QListWidget();
  //  input_files_to_process_list->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);
    input_files_to_process_list->setSortingEnabled(false);

    for (int i = 0; i < list.size(); i++)
    {
        QListWidgetItem* new_item = new QListWidgetItem(list[i]);
        new_item->setTextAlignment(Qt::AlignLeft);
        input_files_to_process_list->addItem(new_item);
    }

    queued_files_layout->addWidget(input_files_to_process_list);

    cout << "New QListWidget Added to Tool Layout" << endl;

    input_files_to_process_lists.push_back(input_files_to_process_list);
    cout << "MainWindow::First name in list: " <<  list.begin()->toStdString() << endl;

    if ((*(list.begin())) != "Brownian"
            && (*(list.begin())) != "LogNormal"
            && (*(list.begin())) != "Gen. Pareto"
            && (*(list.begin())) != "Power Law"
            && (*(list.begin())) != "Exponential"
            && (*(list.begin())) != "CRW"
            && (*(list.begin())) != "LogMCRW"
            && (*(list.begin())) != "GammaMCRW"
            && (*(list.begin())) != "Gamma"
            && (*(list.begin())) != "BootstrapSpeeds"
            && (*(list.begin())) != "BootstrapBoth")
    {

    QStringList parts = list[0].split(QDir::separator());
    dataset_labels.push_back(parts.at(parts.size()-2));
    //dataset_labels[n_datasets] = parts.at(parts.size()-2);
    }
    else
    {
        dataset_labels.push_back(*(list.begin()));
        //dataset_labels[n_datasets] = *(list.begin());
    }

    n_datasets++;

    // setup the efficiency plot legend
/*
if (n_datasets > 0)
{
        time_eff_samples[1]->setName(dataset_labels[0]);
        dist_eff_samples[1]->setName(dataset_labels[0]);
        if (n_datasets > 1)
        {
            time_eff_samples[6]->setName(dataset_labels[1]);
            dist_eff_samples[6]->setName(dataset_labels[1]);

            if (n_datasets > 2)
            {
                time_eff_samples[11]->setName(dataset_labels[2]);
                dist_eff_samples[11]->setName(dataset_labels[2]);
            }
        }
 }
*/
//        timeCustomPlot->replot();
//        distCustomPlot->replot();

    input_files_to_process_list->setMinimumWidth(input_files_to_process_list->sizeHintForColumn(0)/3);

  //  input_files_to_process_list->setMinimumWidth(input_files_to_process_list->sizeHintrColumn(0));
}

void MainWindow::setInputFileProcessed(int v)
{
    for( int i = 0; i < input_files_to_process_lists.size(); i++ )
    {
        if ( v < input_files_to_process_lists[i]->count() )
        {
            input_files_to_process_lists[i]->item(v)->setForeground(Qt::gray);
            break;
        }
        v = v - input_files_to_process_lists[i]->count();
    }
}

void MainWindow::setInputFileProcessing(int v)
{

    cout << "MainWindow::setInputFileProcessing(): v = " << v << endl;
    bool placed = false;

    for( int i = 0; i < input_files_to_process_lists.size(); i++ )
    {
        for( int j = 0; j < input_files_to_process_lists[i]->count(); j++ )
        {
            input_files_to_process_lists[i]->item(j)->setForeground(Qt::black);
        }
    }

    while(!placed)
    {
        int elements_in_this_list = 0;
        int previous_total = 0;

        for( int i = 0; i < input_files_to_process_lists.size(); i++ )
        {
            elements_in_this_list = input_files_to_process_lists[i]->count();
            cout << "MainWindow::setInputFileProcessing(): previous_total = " << previous_total << ", elements in this list: " << elements_in_this_list << endl;

            if ( v < previous_total + elements_in_this_list )
            {
                cout << "MainWindow::setInputFileProcessing(): highlighting element " << v-previous_total << "in list " << i << endl;

                input_files_to_process_lists[i]->setCurrentRow(v-previous_total);
                input_files_to_process_lists[i]->item(v-previous_total)->setSelected(false);
                input_files_to_process_lists[i]->item(v-previous_total)->setForeground(Qt::blue);
                placed = true;
                break;
            }
            previous_total += elements_in_this_list;
        }
        v = v - previous_total;
    }

//    for( int i = 0; i < input_files_to_process_lists.size(); i++ )
//    {
//        if ( v < input_files_to_process_lists[i]->count() )
//        {
//            input_files_to_process_lists[i]->item(v)->setForeground(Qt::blue);
//            input_files_to_process_lists[i]->setCurrentRow(v);
//            break;
//        }
//        v = v - input_files_to_process_lists[i]->count();
//    }
}

void MainWindow::clearInputFilesToProcess()
{
    //input_files_to_process_list->clear();
}


// Format: first value is the sample index (1-5), the second value is the median, (3) lower quartile, (4) upper quartile, (5) min value, (6) max value, remaining values are the outliers.
void MainWindow::setFirstContactSample(FloatVector values)
{
    //cout << "MainWindow::setSample() signal recieved" << endl;

    if (values.size() < 7)
        cout << "MainWindow: too few values in sample: " << values.size() << endl;

    //cout << "MainWindow::plotting at " << values[0] << endl;


    firstContactCustomPlot->setAutoAddPlottableToLegend(false);

int position = values[0];
int setid = values[1];

position -= 1;

int n_samples = first_contact_samples.size();
int difference = position-n_samples+1;
if (difference>0)
    for (int i = 0; i < difference; i++)
    {
        int sample_size = first_contact_samples.size();
        first_contact_samples.push_back(new QCPStatisticalBox(firstContactCustomPlot->xAxis, firstContactCustomPlot->yAxis));
        first_contact_samples[sample_size]->setKey(sample_size);
        first_contact_samples[sample_size]->setOutlierStyle(QCPScatterStyle(QCPScatterStyle::ssPlus, QColor(60, 60, 255, 100), 5));
        first_contact_samples[sample_size]->setMeanStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::black, 8));

        first_contact_samples[sample_size]->setName(dataset_labels[setid]);

        if (setid >= boxBrushes.size())
        {
            QBrush random(QColor(rand()%255,rand()%255,rand()%255,100));
            boxBrushes.push_back(random);
        }

        first_contact_samples[sample_size]->setBrush(boxBrushes[setid]);

        first_contact_x_axis_tick_vector_label << QString::number(sample_size,'f',0);
        first_contact_x_axis_tick_vector << sample_size;
//        cout << "i: " << i << ", position: " << position << ", difference: " << difference << endl;
//        cout << "Added sample: " << sample_size << endl;
        firstContactCustomPlot->addPlottable(first_contact_samples[sample_size]);

        if (isNewFirstContactSetID(setid))
        {
            old_first_contact_setids.push_back(setid);
            first_contact_samples[sample_size]->addToLegend();
        }

    }
//cout << "Updating sample at position: " << position << endl;

first_contact_samples[position]->setMedian(values[2]);
first_contact_samples[position]->setLowerQuartile(values[3]);
first_contact_samples[position]->setUpperQuartile(values[4]);
// not for this histogram -  negative values allowed float min = values[5] < 0? 0: values[5];
first_contact_samples[position]->setMinimum(values[5]);
first_contact_samples[position]->setMaximum(values[6]);
first_contact_samples[position]->setMean(values[7]);

float hopkins = values[8];
float cluster_radius = values[9];
float mean_speed = values[10];
//first_contact_x_axis_tick_vector_label[position] = QString::number(cluster_radius, 'f', 0) + "\n(" + QString::number(hopkins, 'g', 2) + "," + QString::number(mean_speed, 'g', 2) + ")";

//cout << "MainWindow::Cluster Radius: " << x_axis_tick_vector_label[position].toStdString() << endl;

// prepare manual x axis labels:
firstContactCustomPlot->xAxis->setSubTickCount(0);
firstContactCustomPlot->xAxis->setTickLength(0, 4);
firstContactCustomPlot->xAxis->setTickLabelRotation(0);
firstContactCustomPlot->xAxis->setAutoTicks(false);
firstContactCustomPlot->xAxis->setAutoTickLabels(false);
firstContactCustomPlot->xAxis->setTickVector(first_contact_x_axis_tick_vector);
firstContactCustomPlot->xAxis->setTickVectorLabels(first_contact_x_axis_tick_vector_label);


// prepare axes:
firstContactCustomPlot->yAxis->setLabel("Change in Median Time Efficiency (%)"); // ("+time_units_value->text()+")");
firstContactCustomPlot->rescaleAxes();
firstContactCustomPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
//QString s = "Target Cluster Radius ("+spatial_units_value->text()+")\n(Hopkins Statistic, Mean Speed)";
//firstContactCustomPlot->xAxis->setLabel(QString::fromUtf8(s.toStdString().c_str()));

// Prepare the legend
firstContactCustomPlot->legend->setVisible(true);

QFont legendFont = font();  // start out with MainWindow's font..
//legendFont.setPointSize(9); // and make a bit smaller for legend
firstContactCustomPlot->legend->setFont(legendFont);

/*
   time_eff_samples[position]->setMedian(values[2]);
   time_eff_samples[position]->setLowerQuartile(values[3]);
   time_eff_samples[position]->setUpperQuartile(values[4]);
   time_eff_samples[position]->setMinimum(values[5]);
   time_eff_samples[position]->setMaximum(values[6]);
*/
    QVector<double> outliers;

    if (display_outliers)
    for (int i = 11; i < values.size(); i++)
        outliers << values[i];

   first_contact_samples[position]->setOutliers(outliers);

    float max_value = -10000;
    float min_value = 10000;

    /*
    for (int i = 3; i < values.size(); i++)
    {
        if (values[i] > max_value) max_value = values[i];
        if (values[i] < min_value) min_value = values[i];

    }
    */

    max_value = values[6];
    min_value = values[5];

    firstContactCustomPlot->yAxis->setRange(min_value - (abs(min_value-max_value))*0.1, max_value + (abs(min_value-max_value)*0.1));

    firstContactCustomPlot->rescaleAxes();
    firstContactCustomPlot->yAxis->scaleRange(1.1, firstContactCustomPlot->yAxis->range().center());
    firstContactCustomPlot->xAxis->scaleRange(1.1, firstContactCustomPlot->xAxis->range().center());

    firstContactCustomPlot->replot();
}

// Format: first value is the sample index (1-5), the second value is the median, (3) lower quartile, (4) upper quartile, (5) min value, (6) max value, remaining values are the outliers.
void MainWindow::setTimeEfficiencySample(FloatVector values)
{
    //cout << "MainWindow::setSample() signal recieved" << endl;

    if (values.size() < 7)
        cout << "MainWindow: too few values in sample: " << values.size() << endl;

    //cout << "MainWindow::plotting at " << values[0] << endl;


    timeCustomPlot->setAutoAddPlottableToLegend(false);

int position = values[0];
int setid = values[1];

position -= 1;

int n_samples = time_eff_samples.size();
int difference = position-n_samples+1;
if (difference>0)
    for (int i = 0; i < difference; i++)
    {
        int sample_size = time_eff_samples.size();
        time_eff_samples.push_back(new QCPStatisticalBox(timeCustomPlot->xAxis, timeCustomPlot->yAxis));
        time_eff_samples[sample_size]->setKey(sample_size);
        time_eff_samples[sample_size]->setOutlierStyle(QCPScatterStyle(QCPScatterStyle::ssPlus, QColor(60, 60, 255, 100), 5));
        time_eff_samples[sample_size]->setMeanStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::black, 8));

        time_eff_samples[sample_size]->setName(dataset_labels[setid]);

        if (setid >= boxBrushes.size())
        {
            QBrush random(QColor(rand()%255,rand()%255,rand()%255,100));
            boxBrushes.push_back(random);
        }

        time_eff_samples[sample_size]->setBrush(boxBrushes[setid]);

        time_x_axis_tick_vector_label << QString::number(sample_size,'f',0);
        time_x_axis_tick_vector << sample_size;
//        cout << "i: " << i << ", position: " << position << ", difference: " << difference << endl;
//        cout << "Added sample: " << sample_size << endl;
        timeCustomPlot->addPlottable(time_eff_samples[sample_size]);

        if (isNewTimeSetID(setid))
        {
            old_timesetids.push_back(setid);
            time_eff_samples[sample_size]->addToLegend();
        }

    }
//cout << "Updating sample at position: " << position << endl;

time_eff_samples[position]->setMedian(values[2]);
time_eff_samples[position]->setLowerQuartile(values[3]);
time_eff_samples[position]->setUpperQuartile(values[4]);
float min = values[5] < 0? 0: values[5];
time_eff_samples[position]->setMinimum(min);
time_eff_samples[position]->setMaximum(values[6]);
time_eff_samples[position]->setMean(values[7]);

float hopkins = values[8];
float cluster_radius = values[9];
float time_used = values[10];
time_x_axis_tick_vector_label[position] = QString::number(cluster_radius, 'f', 0) + "\n(" + QString::number(hopkins, 'g', 2) + "," + QString::number(time_used, 'g', 2) + ")";

//cout << "MainWindow::Cluster Radius: " << x_axis_tick_vector_label[position].toStdString() << endl;

// prepare manual x axis labels:
timeCustomPlot->xAxis->setSubTickCount(0);
timeCustomPlot->xAxis->setTickLength(0, 4);
timeCustomPlot->xAxis->setTickLabelRotation(0);
timeCustomPlot->xAxis->setAutoTicks(false);
timeCustomPlot->xAxis->setAutoTickLabels(false);
timeCustomPlot->xAxis->setTickVector(time_x_axis_tick_vector);
timeCustomPlot->xAxis->setTickVectorLabels(time_x_axis_tick_vector_label);

// prepare axes:
timeCustomPlot->yAxis->setLabel("Time Efficiency (Targets/"+time_units_value->text()+")");
timeCustomPlot->rescaleAxes();
timeCustomPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
QString s = "Target Cluster Radius ("+spatial_units_value->text()+")\n(Hopkins Statistic, Mean Time Used)";
timeCustomPlot->xAxis->setLabel(QString::fromUtf8(s.toStdString().c_str()));

// Prepare the legend
timeCustomPlot->legend->setVisible(true);

QFont legendFont = font();  // start out with MainWindow's font..
//legendFont.setPointSize(9); // and make a bit smaller for legend
timeCustomPlot->legend->setFont(legendFont);

/*
   time_eff_samples[position]->setMedian(values[2]);
   time_eff_samples[position]->setLowerQuartile(values[3]);
   time_eff_samples[position]->setUpperQuartile(values[4]);
   time_eff_samples[position]->setMinimum(values[5]);
   time_eff_samples[position]->setMaximum(values[6]);
*/
    QVector<double> outliers;

    if (display_outliers)
    for (int i = 11; i < values.size(); i++)
        outliers << values[i];

   time_eff_samples[position]->setOutliers(outliers);

    float max_value = -10000;
    float min_value = 10000;

    /*
    for (int i = 3; i < values.size(); i++)
    {
        if (values[i] > max_value) max_value = values[i];
        if (values[i] < min_value) min_value = values[i];

    }
    */

    max_value = values[6];
    min_value = values[5];

    timeCustomPlot->yAxis->setRange(min_value - (abs(min_value-max_value))*0.1, max_value + (abs(min_value-max_value)*0.1));

    timeCustomPlot->rescaleAxes();
    timeCustomPlot->yAxis->scaleRange(1.1, timeCustomPlot->yAxis->range().center());
    timeCustomPlot->xAxis->scaleRange(1.1, timeCustomPlot->xAxis->range().center());

    timeCustomPlot->replot();
}

void MainWindow::setDistanceEfficiencySample(FloatVector values)
{
    //cout << "MainWindow::setSample() signal recieved" << endl;

    if (values.size() < 7)
        cout << "MainWindow: too few values in sample: " << values.size() << endl;

    //cout << "MainWindow::plotting at " << values[0] << endl;

    distCustomPlot->setAutoAddPlottableToLegend(false);

int position = values[0];
int setid = values[1];

position -= 1;

int n_samples = dist_eff_samples.size();
int difference = position-n_samples+1;
if (difference>0)
    for (int i = 0; i < difference; i++)
    {
        int sample_size = dist_eff_samples.size();
        dist_eff_samples.push_back(new QCPStatisticalBox(distCustomPlot->xAxis, distCustomPlot->yAxis));
        dist_eff_samples[sample_size]->setKey(sample_size);
        dist_eff_samples[sample_size]->setOutlierStyle(QCPScatterStyle(QCPScatterStyle::ssPlus, QColor(60, 60, 255, 100), 5));
        dist_eff_samples[sample_size]->setMeanStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::black, 8));

        dist_eff_samples[sample_size]->setName(dataset_labels[setid]);

        if (setid >= boxBrushes.size())
        {
            QBrush random(QColor(rand()%255,rand()%255,rand()%255,100));
            boxBrushes.push_back(random);
        }

        dist_eff_samples[sample_size]->setBrush(boxBrushes[setid]);

        dist_x_axis_tick_vector_label << QString::number(sample_size);
        dist_x_axis_tick_vector << sample_size;
        cout << "i: " << i << ", position: " << position << ", difference: " << difference << endl;
        cout << "Added sample: " << sample_size << endl;
        distCustomPlot->addPlottable(dist_eff_samples[sample_size]);

        if (isNewDistSetID(setid))
        {
            old_distsetids.push_back(setid);
            dist_eff_samples[sample_size]->addToLegend();
        }

    }
//cout << "Updating sample at position: " << position << endl;

dist_eff_samples[position]->setMedian(values[2]);
dist_eff_samples[position]->setLowerQuartile(values[3]);
dist_eff_samples[position]->setUpperQuartile(values[4]);
float min = values[5] < 0? 0: values[5];
dist_eff_samples[position]->setMinimum(min);
dist_eff_samples[position]->setMaximum(values[6]);
dist_eff_samples[position]->setMean(values[7]);

float hopkins = values[8];
float cluster_radius = values[9];
float dist_used = values[10];
dist_x_axis_tick_vector_label[position] = QString::number(cluster_radius, 'f', 0) + "\n(" + QString::number(hopkins, 'g', 2) + "," + QString::number(dist_used, 'g', 2) + ")";

// prepare manual x axis labels:
distCustomPlot->xAxis->setSubTickCount(0);
distCustomPlot->xAxis->setTickLength(0, 4);
distCustomPlot->xAxis->setTickLabelRotation(0);
distCustomPlot->xAxis->setAutoTicks(false);
distCustomPlot->xAxis->setAutoTickLabels(false);
distCustomPlot->xAxis->setTickVector(dist_x_axis_tick_vector);
distCustomPlot->xAxis->setTickVectorLabels(dist_x_axis_tick_vector_label);


// prepare axes:
distCustomPlot->yAxis->setLabel("Distance Efficiency (Targets/"+spatial_units_value->text()+")");
distCustomPlot->rescaleAxes();
distCustomPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
QString s = "Target Cluster Radius ("+spatial_units_value->text()+")\n(Hopkins Statistic, Mean Dist. Used)";
distCustomPlot->xAxis->setLabel(QString::fromUtf8(s.toStdString().c_str()));

// Prepare the legend
distCustomPlot->legend->setVisible(true);

QFont legendFont = font();  // start out with MainWindow's font..
//legendFont.setPointSize(9); // and make a bit smaller for legend
distCustomPlot->legend->setFont(legendFont);

/*
   time_eff_samples[position]->setMedian(values[2]);
   time_eff_samples[position]->setLowerQuartile(values[3]);
   time_eff_samples[position]->setUpperQuartile(values[4]);
   time_eff_samples[position]->setMinimum(values[5]);
   time_eff_samples[position]->setMaximum(values[6]);
*/
    QVector<double> outliers;

    if (display_outliers)
    for (int i = 11; i < values.size(); i++)
        outliers << values[i];

   dist_eff_samples[position]->setOutliers(outliers);

    float max_value = -10000;
    float min_value = 10000;

    /*
    for (int i = 3; i < values.size(); i++)
    {
        if (values[i] > max_value) max_value = values[i];
        if (values[i] < min_value) min_value = values[i];

    }
    */

    max_value = values[6];
    min_value = values[5];

    distCustomPlot->yAxis->setRange(min_value - (abs(min_value-max_value))*0.1, max_value + (abs(min_value-max_value)*0.1));

    distCustomPlot->rescaleAxes();
    distCustomPlot->yAxis->scaleRange(1.1, distCustomPlot->yAxis->range().center());
    distCustomPlot->xAxis->scaleRange(1.1, distCustomPlot->xAxis->range().center());

    distCustomPlot->replot();
}


//// Format: first value is the sample index (1-5), the second value is the median, (3) lower quartile, (4) upper quartile, (5) min value, (6) max value, remaining values are the outliers.
//void MainWindow::setDistanceEfficiencySample(FloatVector values)
//{
//    cout << "MainWindow::setSample() signal recieved" << endl;

//    if (values.size() < 7)
//        cout << "MainWindow: too few values in sample: " << values.size() << endl;

//    //cout << "MainWindow::plotting at " << values[0] << endl;

// int position = values[0];

// position -= 1;

// int n_samples = dist_eff_samples.size();
// int difference = position-n_samples+1;
// if (difference>0)
// for (int i = 0; i < difference; i++)
// {
//     dist_eff_samples.push_back(new QCPStatisticalBox(distCustomPlot->xAxis, distCustomPlot->yAxis));
//     x_axis_tick_vector_label.push_back("New");
// }

// dist_eff_samples[position]->setMedian(values[2]);
// dist_eff_samples[position]->setLowerQuartile(values[3]);
// dist_eff_samples[position]->setUpperQuartile(values[4]);
// dist_eff_samples[position]->setMinimum(values[5]);
// dist_eff_samples[position]->setMaximum(values[6]);
// dist_eff_samples[position]->setMean(values[7]);
// float hopkins = values[8];
// float cluster_radius = values[9];

// /*
//  dist_eff_samples[position]->setMedian(values[2]);
//  dist_eff_samples[position]->setLowerQuartile(values[3]);
//  dist_eff_samples[position]->setUpperQuartile(values[4]);
//  dist_eff_samples[position]->setMinimum(values[5]);
//  dist_eff_samples[position]->setMaximum(values[6]);
//*/
//     x_axis_tick_vector_label[position] = QString::number(cluster_radius, 'f', 0) + "\n(" + QString::number(hopkins, 'g', 2) + ")";


// distCustomPlot->xAxis->setTickVectorLabels(x_axis_tick_vector_label);


//    QVector<double> outliers;

//    if (display_outliers)
//    for (int i = 10; i < values.size(); i++)
//        outliers << values[i];

//  dist_eff_samples[position]->setOutliers(outliers);

//    float max_value = -10000;
//    float min_value = 10000;
//    /*
//    for (int i = 3; i < values.size(); i++)
//    {
//        if (values[i] > max_value) max_value = values[i];
//        if (values[i] < min_value) min_value = values[i];

//    }
//    */

//    max_value = values[6];
//    min_value = values[5];

//    distCustomPlot->yAxis->setRange(min_value - (abs(min_value-max_value))*0.1, max_value + (abs(min_value-max_value)*0.1));

//    distCustomPlot->rescaleAxes();

//    distCustomPlot->replot();

//    //distCustomPlot->rescaleValueAxis();

//}

void MainWindow::save()
{
    cout << "Save called." << endl;

    QString dir = QFileDialog::getExistingDirectory(this, tr("Select Directory"),
                                                "~",
                                                QFileDialog::ShowDirsOnly
                                                | QFileDialog::DontResolveSymlinks);

    setOutputDirPath(dir);


    //QString path = QFileDialog::getSaveFileName(this, tr("Set Output Directory"), default_save_path, tr("Portable Document Format File (*.pdf)"));
}

void MainWindow::updateTime()
{
    int elapsed_milliseconds = elapsed_time.elapsed();

    int secs = elapsed_milliseconds / 1000;
    int mins = (secs / 60) % 60;
    int hours = (secs / 3600);
    secs = secs % 60;
    elapsed_time_value->setText(QString("%1:%2:%3")
    .arg(hours, 2, 10, QLatin1Char('0'))
    .arg(mins, 2, 10, QLatin1Char('0'))
    .arg(secs, 2, 10, QLatin1Char('0')) );

    int remaining_secs = est_remaining_time_secs;
    //cout << "MainWindow::updateTime(): " << est_remaining_time_secs << endl;
    int remaining_mins = (remaining_secs / 60) % 60;
    int remaining_hours = (remaining_secs / 3600);
    remaining_secs = remaining_secs % 60;

    est_remaining_time_value->setText(QString("%1:%2:%3")
    .arg(remaining_hours, 2, 10, QLatin1Char('0'))
    .arg(remaining_mins, 2, 10, QLatin1Char('0'))
    .arg(remaining_secs, 2, 10, QLatin1Char('0')) );
}

void MainWindow::setTimeUnits(QString v)
{
    timeCustomPlot->yAxis->setLabel("Time Efficiency (Targets/"+time_units_value->text()+")");
    time_efficiency_units_label->setText(time_units_value->text());
    setTimeResolution(time_resolution);
    setTimeExpended(dist_expended);
    timeCustomPlot->replot();
    distCustomPlot->replot();
    speed_histogram->setXUnits(spatial_units_value->text() + "/" + time_units_value->text());
    speed_histogram->repaint();
}

void MainWindow::setSpatialUnits(QString v)
{
    distCustomPlot->yAxis->setLabel("Distance Efficiency (Targets/"+spatial_units_value->text()+")");
    speed_histogram->setXUnits(spatial_units_value->text() + "/" + time_units_value->text());
    speed_histogram->repaint();
    dist_efficiency_units_label->setText(spatial_units_value->text());
    setDensity(density);
    setDistExpended(dist_expended);
    QString s = "Target Cluster Radius ("+spatial_units_value->text()+")\n(Hopkins Statistic)";
    timeCustomPlot->xAxis->setLabel(s);
    distCustomPlot->xAxis->setLabel(s);
    timeCustomPlot->replot();
    distCustomPlot->replot();
}

void MainWindow::addDataset()
{
    cout << "MainWindow: addDataset() called" << endl;
    if (getModel()->getSearchType() == 0)
    {
        OpenFile();
    }
    else
    {
        QStringList idealized;

        switch(getModel()->getSearchType())
        {
        case 1: idealized.push_back("Brownian"); break;
        case 2: idealized.push_back("LogNormal"); break;
        case 3: idealized.push_back("Gen. Pareto"); break;
        case 4: idealized.push_back("Power Law"); break;
        case 5: idealized.push_back("Exponential"); break;
        case 6: idealized.push_back("Gamma"); break;
        case 7: idealized.push_back("CRW"); break;
        case 8: idealized.push_back("LogMCRW"); break;
        case 9: idealized.push_back("GammaMCRW"); break;
        case 10: idealized.push_back("BootstrapSpeeds"); break;
        case 11: idealized.push_back("BootstrapBoth"); break;
        }

        getModel()->queueInputFiles( idealized );
    }
}

QStringList MainWindow::OpenFile()
{
    QString default_path = QString::fromStdString(getModel()->getWorkingPath());

#if defined(__APPLE__)
    default_path += "/MyFile.txt";
#endif

    QStringList fileNames = QFileDialog::getOpenFileNames(this, tr("Open File"), default_path, tr("Comma Separated Values File (*.csv)"));

    if (!fileNames.isEmpty())
        getModel()->queueInputFiles( fileNames );

    return fileNames;
}

Model* MainWindow::getModel()
{
    return vis->getModel();
}

void MainWindow::setDisplayOutliers(bool v)
{
    display_outliers = v;
    distCustomPlot->rescaleAxes();
    timeCustomPlot->rescaleAxes();
    distCustomPlot->replot();
    timeCustomPlot->replot();
}

bool MainWindow::isNewTimeSetID(int setid)
{
    for (int i = 0; i < old_timesetids.size(); i++)
    {
        if (old_timesetids[i] == setid) return false;
    }

    return true;
}

bool MainWindow::isNewFirstContactSetID(int setid)
{
    for (int i = 0; i < old_first_contact_setids.size(); i++)
    {
        if (old_first_contact_setids[i] == setid) return false;
    }

    return true;
}

bool MainWindow::isNewDistSetID(int setid)
{
    for (int i = 0; i < old_distsetids.size(); i++)
    {
        if (old_distsetids[i] == setid) return false;
    }

    return true;
}

void MainWindow::setDistScaleType(bool v)
{
    cout << "MainWindow::setDistScaleType() called with: " << v << endl;
    QCPAxis::ScaleType type;
    if (v)
        type = QCPAxis::stLogarithmic;
    else
        type = QCPAxis::stLinear;

    distCustomPlot->yAxis->setScaleType(type);
    distCustomPlot->rescaleAxes();
    distCustomPlot->replot();

}

void MainWindow::setTimeScaleType(bool v)
{
    cout << "MainWindow::setTimeScaleType() called with: " << v << endl;
    QCPAxis::ScaleType type;
    if (v)
        type = QCPAxis::stLogarithmic;
    else
        type = QCPAxis::stLinear;

    timeCustomPlot->yAxis->setScaleType(type);
    timeCustomPlot->rescaleAxes();
    timeCustomPlot->replot();

}

void MainWindow::setOutputDirPath(QString path)
{
    default_save_path = path;
    output_dir_path_value->setText(path);
    emit setModelSavePath(path);
}

void MainWindow::modelFinished()
{
    stringstream time_ss, dist_ss, diff_ss;
    time_ss << default_save_path.toStdString() << "/time.pdf";
    dist_ss << default_save_path.toStdString() << "/dist.pdf";
    diff_ss << default_save_path.toStdString() << "/diff.pdf";

    QString time_path = QString::fromStdString(time_ss.str());
    QString dist_path = QString::fromStdString(dist_ss.str());
    QString diff_path = QString::fromStdString(diff_ss.str());

    timeCustomPlot->savePdf(time_path);
    distCustomPlot->savePdf(dist_path);
    firstContactCustomPlot->savePdf(diff_path);

    cout << "MainWindow: saved time efficiency plot to " << time_path.toStdString() << endl;
    cout << "MainWindow: saved distance efficiency plot to " << dist_path.toStdString() << endl;
    cout << "MainWindow: saved difference efficiency plot to " << diff_path.toStdString() << endl;
}
