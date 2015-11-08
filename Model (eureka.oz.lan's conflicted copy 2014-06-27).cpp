//Weighted average by distance covered
//report the volume.

#include "MainWindow.h"
#include "QDebugStream.h"
#include "Model.h"
#include "Sphere.h"
#include "RectangularPrism.h"
#include "Searcher.h"
#include "Target.h"
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <QLinkedList>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <QThread>
#include <numeric>
#include "Worker.h"
#include <QVector>
#include <QFileInfo>
#include <algorithm>
#include "ALGLIB/statistics.h"

using namespace std;
using boost::lexical_cast;
using namespace boost::math;
using namespace alglib;

Model::Model( float spr, float cr, int nc, int ntc, int sp, bool t, string app_path )
{
  // T cell simulation paramters
  // Search Sphere size 1mm^3
  // Destructive: true
  // Simulation Duration: 1hr, 4hr, 12hr
  // Cellular density:
  // 50% e filled with DCs and Tcells
  // DCs 2200 micrometer volume
  // 2-5% of Cellular density is made up of DCs
  // . Miller, M.J., S.H. Wei, M.D. Cahalan, and I. Parker. 2004.
  // T cell repertoire scanning is promoted by dynamic dendritic
  // cell behavior and random T cell motility in the lymph node.
  // Proc. Natl. Acad. Sci. USA. 101:998–1003
  //
  // 1mm^3 = 10^9um^3
  // DC has volume 1000um^3
  // 10^6um^3
  //
  // 1mm = 10^(-3)um
  // 1mm^3 = 10^(-3)^3um^3 = 10^(-9)um^3
  //
  // Diameter = 10um
  //
  // 100 DCs per edge.
  //
  // 100^3 DCs fill the LNs volume.
  //
  // 10^6 DCs fill the LNs volume.
  //
  // 0.5% of 10^6 = 5*10^(-3) (Percentage of cells that are DCs)
  //
  // 5*10^(-3)*10^6 = 5*10^3
  //
  // 5000 DCs in the Lymph Node.
  //
  // 70 groups of 70

  // Each DC has radius 0.01 (10^(-2)) of the sphere.


    reference_volume = 300*300*70; // micrometers
    scale_number_of_targets = 1.0;

    gui_ready = false;

    targets_mutex = false;
    targets_mutex_seen = false;

    searchers_mutex = false;
    searchers_mutex_seen = false;

    paths_mutex = false;
    paths_mutex_seen = false;
    num_targets_found = 0;
    stopped = true;

    this->time_resolution = 1;

    iTimer = new QTimer(this);
    iTimer->setInterval(1000/30);//75);
    QObject::connect(iTimer, SIGNAL(timeout()), this, SIGNAL(updateNeeded()));
    iTimer->start();

  targets_with_replacement = false;
  time_efficiency = t;
  max_targets = 1000;
  search_pattern = sp;
  detection_radius = spr;
  cluster_radius = cr;
  n_clusters = nc;
  targets_per_cluster = ntc;

  time_used = 0;
  distance_used = 0;
  
  time = 0;
    max_time = 0;
    max_dist = 0;
    n_targets_found = 0;

    float radius = 1.0;

    float x_min = 0;
    float x_max = 100;

    float y_min = 0;
    float y_max = 100;

    float z_min = 0;
    float z_max = 100;

    readSettings();

    space = new RectangularPrism(x_min, x_max, y_min, y_max, z_min, z_max, n_clusters*targets_per_cluster);//new Sphere( radius, n_clusters*targets_per_cluster, detection_radius, search_pattern );

    // Create searchers
    /*
    int n_searchers = space->getNumSearchers();
    Searcher* searchers = space->getSearchers();

    //   float x = 0, y = 0, z = 0;

    //for ( int i = 0; i < n_searchers; i++ )
    //  {
    //    searchers[i].setXPos(x);
    //    searchers[i].setYPos(y);
    //    searchers[i].setZPos(z);
    //  }

        float x, y, z;

        int i = 0;
        while( i < n_searchers )
        {
            x = (rand()%2000/1000.0-1.0f);
            y = (rand()%2000/1000.0-1.0f);
            z = (rand()%2000/1000.0-1.0f);

            if ( x*x + y*y + z*z < 1 )
            {
                searchers[i].setXPos(x);
                searchers[i].setYPos(y);
                searchers[i].setZPos(z);
                searchers[i].setPsiAngle(0);
                searchers[i].setThetaAngle(0);
                searchers[i].setPhiAngle(0);
                searchers[i].setRadius(0.005);
                i++;
            }
        }
*/

    //placeTargets();

}

void Model::setClusterRadius(float v)
{
    //cout << "Model: cluster radius set to " << v << endl;

    if (v < 5)
    {
        cout << "Model: cluster radius too small" << endl;
        return;
    }

    cluster_radius = v;

    placeTargets();

    //emit updateNeeded();

}

SearchSpace* Model::getSearchSpace()
{
    return space;
}

float Model::getClusterRadius()
{
    return cluster_radius;
}

void Model::generate()
{

    //cout << "Max Time: " << max_time << endl;

    //cout << "Generate called in thread " << pthread_self() << endl;

    float time = 0;

    vector<float> velocities;

    distance_used = 0;
    time_used = 0;

    while ( time < max_time )
    {
        //qApp->processEvents();

        ////cout << "Generated velocity" << endl;
        ////cout << "Updating in thread " << pthread_self() << ". Step " << time << endl;


        vector<float> new_velocities = space->updateSearcherPositions().toStdVector();


        velocities.insert(velocities.end(), new_velocities.begin(), new_velocities.end());

        time = velocities.size()*time_resolution;

        //calcTargetsFound();
        //emit updateNeeded();
    }

    time_used = space->getTotalTimeTravelled();
    distance_used = space->getTotalDistanceTravelled();

    // Calculate the histogram values
        // Get the max and min values in velocities. max - min is the velocity range
        // Create 10 bins with width bin range 1/n_bins of the range.
        // Loop over all the velocities. If the value falls in the bin range increment histogram

    QVector<int> histogram;

    int n_bins = 100;
    float bin_mins[n_bins];

    // initialize the histogram vector
    for (int i = 0; i < n_bins; i++)
    {
            histogram.append(0);
        }

    if (velocities.empty())
    {
        cout << "Model::generateSearch(): no velocities found." << endl;
        return;
    }

    // Get velocity min and max
    vector<float>::const_iterator it;

    it = max_element(velocities.begin(), velocities.end());
    float max_velocity = *it;

    it = min_element(velocities.begin(), velocities.end());

    float min_velocity = *it;

    float range = max_velocity - min_velocity;

    float bin_size = range/n_bins;

    // initialize the bin_mins array
    float cumulation = 0;
    for (int i = 0; i < n_bins; i++)
    {
        bin_mins[i] = cumulation;
        cumulation += bin_size;
    }

    for(std::vector<float>::iterator it = velocities.begin(); it != velocities.end(); ++it)
    {
        float v = *it;

        for (int i = 0; i < n_bins; i++)
        {
            if (v < bin_mins[i])
            {
                histogram[i]++;
                break;
            }
        }

    }

    //cout << "Number of velocities: " << velocities.size() <<  endl;
    //cout << "Bin Min Values: ";
    //for ( int i = 0; i < n_bins; i++) cout << bin_mins[i] << " ";
    //cout << endl;

    //cout << "Bin Counts: ";
    //for ( int i = 0; i < histogram.size(); i++) cout << histogram[i] << " ";
    //cout << endl;


    emit updateTimeExpended(time_used);
    emit updateDistExpended(distance_used);
    emit updateHistogramValues(histogram);
    emit updateHistogramMaxValue(cumulation);
    qApp->processEvents();

}

void Model::run()
{
    //cout << "Run called" << endl;
    //iTimer->start();
    //sleep(250);
    generate();
   // emit updateNeeded();

    /*
    if (stopped)
    {

        return;
    }

    if (mode == GENERATE)
    {
    //while(1)
    //{

    time++;
   // iTimer->start();

    //if (time == 100000) exit(1);

    //Searcher* searchers = space->getSearchers();

    space->updateSearcherPositions( 1 );

    //Searcher* searchers = space->getSearchers();
    //for (int i = 0; i < space->getNumSearchers(); i++)
    //  cout << time << "," << searchers[i].getXPos() << "," << searchers[i].getYPos() << "," << searchers[i].getZPos() << endl;

    int n_targets_found = 0;
    Target* targets = space->getTargets();
    // Check whether searchers have found targets
      for ( int j = 0; j < space->getNumTargets(); j++ )
	  {
	   
	    if (targets[j].isFound())
	      n_targets_found++;
	  }

      num_targets_found = n_targets_found;

      cout << "Targets Found: " << num_targets_found << endl;

    //cout << time << ": " << n_targets_found << endl;

    }
    else if (mode == RETRACE)
    {
        // Calculate target interactions
        int n_searchers = space->getNumSearchers();
        Searcher* searchers = space->getSearchers();

        for ( int i = 0; i < n_searchers; i++ )
        {
            Searcher s = searchers[i];
            vector<Coordinate*> path = s.getPath();

            for (vector<Coordinate*>::iterator it = path.begin(); it != path.end(); ++it)
            {
                    Coordinate* start = *it;
                    Coordinate* end = *it;
                    space->testpath( start->getX(), start->getY(), start->getZ(), end->getX(), end->getY(), end->getZ() );
            }
        }
    }


    iTimer->start();
    //}
    */
}

void Model::stop()
{
    stopped = true;
}

void Model::reset()
{
    removeTargets();
    placeTargets();
}

void Model::start()
{
    stopped = false;
}

Model::~Model()
{
    //if (worker) worker->abort();
    //if (workerThread) workerThread->wait();
    //if (workerThread) delete workerThread;
    //if (worker) delete worker;
    //delete space;
    saveSettings();
    //cout << "Model: destructor called" << endl;
}


// This function reads a file and loads each track into a searcher. The searchers then follow the paths.
// Should be able to save tracks to csv from within PowerSearch
void Model::addSearchers(QString fp)
{

    int n_filtered = 0;

    std::string filepath = fp.toUtf8().constData();

    //cout << "Reading file " << filepath << endl;

    QFileInfo info(QString::fromStdString(filepath));

    working_path = info.dir().absolutePath().toStdString();

    cout << "Model: working_path set to " << working_path << endl;

    vector<Searcher*> searcher_v;
    Searcher* current_searcher = NULL;
    int id;
    float x, y, z, t, max_x = -10000, max_y = -10000, max_z = -10000, min_x = 10000, min_y = 10000, min_z = 10000, max_t = -10000, min_t = 10000;
    ifstream infile(filepath.c_str());

    if (!infile) {
      cerr << "Can't open input file." << endl;
    }
    else
    {
        //cout << "File opened successfully." << endl;
    }

    string line;

        //while( infile >> id >> x >> y >> z >> t )

        int n_data_points_for_this_id = 0;

        int current_line_number = 0;

        int path_length_threshold = 3;

        while( infile >> line )
        {
            current_line_number++;
            qApp->processEvents();
            //cout << "Read: " << line << endl;
            std::vector<std::string> strs;
            boost::split(strs, line, boost::is_any_of(", "));

            if ( strs.size() != 5 )
            {

                parseError("Parse Error: need at least five valuse per line. In " + fp + ". Line " + QString::number(current_line_number));
                return;
            }

            id = lexical_cast<int>(strs.at(0));
            x = lexical_cast<float>(strs.at(1));
            y = lexical_cast<float>(strs.at(2));
            z = lexical_cast<float>(strs.at(3));
            t = lexical_cast<float>(strs.at(4));

            //z = abs(z);

            if (x > max_x) max_x = x;
            if (y > max_y) max_y = y;
            if (z > max_z) max_z = z;
            if (t > max_t) max_t = t;

            if (x < min_x) min_x = x;
            if (y < min_y) min_y = y;
            if (z < min_z) min_z = z;
            if (t < min_t) min_t = t;

            n_data_points_for_this_id++;

                if ( current_searcher == NULL || current_searcher->getID() != id )
                {
                    //cout << "Data points" << n_data_points_for_this_id << endl;

                    if (current_searcher != NULL)
                    {
                        if (current_searcher->getPath().empty())
                        {
                            cout << "Too little data for track id" << endl;
                            parseError("Parse Error: need at least two poitions per track id. In " + fp + ". Line " + QString::number(current_line_number));
                            return;
                        }
                        if (current_searcher->getPath().size() > path_length_threshold)
                        {
                            // Filter by min path length
                            searcher_v.push_back( current_searcher );
                        }
                        else
                        {
                            cout << "Model::filtered searcher with path length of less than " << path_length_threshold << ". Filterd: " << ++n_filtered << endl;
                        }

                    }

                    n_data_points_for_this_id = 0;
                    // cout << "Creating new searcher with ID " << id << endl;
                    current_searcher = new Searcher();
                    current_searcher->setID(id);


                        current_searcher->setStartingXPos(x);
                        current_searcher->setStartingYPos(y);
                        current_searcher->setStartingZPos(z);
                        current_searcher->setStartingTime(t);
                }



            // cout << "Read < " << id << " " << x << " " << y << " " << z << " " << t << " >" << endl;
            Coordinate* c = new Coordinate(x, y, z, t);
            current_searcher->addToPath(c);
        }


//cout << "Finished reading data file." << endl;

    infile.close();



  //  reset();

    int n_searchers = searcher_v.size();

    cout << "Read " << n_searchers << " searchers." << endl;

    Searcher* searchers = new Searcher[n_searchers];

    int i = 0;
    for (std::vector<Searcher*>::iterator it = searcher_v.begin() ; it != searcher_v.end(); ++it)
    {
        searchers[i].setID((*it)->getID());
        vector<Coordinate*> p = (*it)->getPath();
        if (p.size() >= 2) searchers[i].setPath(p);
        else cout << "Model: searcher path too short" << endl;
        i++;

    }

    //cout << "x_max=" << max_x << " y_max=" << max_y << " z_max=" << max_z << endl;

    time_resolution = min_t;
    space->setTimeResolution(time_resolution);
    emit updateTimeResolution(time_resolution);
    qApp->processEvents();

    // Scale paths
   // for (int i = 0; i < n_searchers; i++)
   // {
   //     searchers[i].scalePath(min_x, max_x, min_y, max_y, min_z, max_z);
   // }


    // Set initial position to the first position read
    for (int i = 0; i < n_searchers; i++)
    {
        vector<Coordinate*> path = searchers[i].getPath();
        if (!path.empty())
        {
            Coordinate* start_position = path.at(0);
            searchers[i].setXPos( start_position->getX() );
            searchers[i].setYPos( start_position->getY() );
            searchers[i].setZPos( start_position->getZ() );
        }
        else
        {
           emit parseError("Parse Error: Searcher found with empty path.");
           return;
        }
    }

    removeSearchers(); // GUI thread safe
    getSearchSpace()->setSearchers( searchers, n_searchers );
    //cout << "Added " << getSearchSpace()->getNumSearchers() << " Searchers" << endl;

    getSearchSpace()->setBounds( min_x, min_y, min_z, max_x, max_y, max_z );
    emit updateVolume(getSearchSpace()->getVolume());

    scale_number_of_targets = getSearchSpace()->getVolume()/reference_volume;


    // Calculate total time
    // Calculate total distance

    float total_dist = 0;
    float total_time = 0;

    float velocity = 0;
    // Calculate the histogram values
        // Get the max and min values in velocities. max - min is the velocity range
        // Create 10 bins with width bin range 1/n_bins of the range.
        // Loop over all the velocities. If the value falls in the bin range increment histogram

    QVector<int> histogram;

    int n_bins = 100;
    float bin_mins[n_bins];

    // initialize the histogram vector
    for (int i = 0; i < n_bins; i++)
    {
            histogram.append(0);
        }


    vector<float> velocities;

    for ( int i = 0; i < n_searchers; i++ )
    {
        Searcher s = searchers[i];
        vector<Coordinate*> path = s.getPath();

        for (int j = 0; j < path.size()-1; j++)
        {

                //cout << "j = " << j << " of " << path.size() << endl;
                Coordinate* start = path[j];
                Coordinate* end = path[j+1];

                float velocity = space->norm(start->getX()-end->getX(),start->getY()-end->getY(),start->getZ()-end->getZ());

                velocities.push_back(velocity);

                total_dist += velocity;
        }
    }

    // Get velocity min and max
    vector<float>::const_iterator it;
    it = max_element(velocities.begin(), velocities.end());
    float max_velocity = *it;

    it = min_element(velocities.begin(), velocities.end());
    float min_velocity = *it;

    float range = max_velocity - min_velocity;

    float bin_size = range/n_bins;

    // initialize the bin_mins array
    float cumulation = 0;
    for (int i = 0; i < n_bins; i++)
    {
        bin_mins[i] = cumulation;
        cumulation += bin_size;
    }

    for(std::vector<float>::iterator it = velocities.begin(); it != velocities.end(); ++it)
    {
        float v = *it;

        for (int i = 0; i < n_bins; i++)
        {
            if (v < bin_mins[i])
            {
                histogram[i]++;
                break;
            }
        }

    }

    //cout << "Time resolution: " << time_resolution << endl;
    //cout << "Number of velocities: " << velocities.size() <<  endl;
    //cout << "Bin Min Values: ";
    //for ( int i = 0; i < n_bins; i++) cout << bin_mins[i] << " ";
    //cout << endl;

    total_time = velocities.size()*time_resolution;

    //cout << "Bin Counts: ";
    //for ( int i = 0; i < histogram.size(); i++) cout << histogram[i] << " ";
    //cout << endl;


    max_time = total_time; // max_time is used when tracks are generated to bound how long they search, ie they search for as long as the real T cells were observed to search
    max_dist = total_dist; // same for max_dist

    //cout << "Total time observed: " << max_time << " seconds." << endl;

    placeTargets();
    calcTargetsFound();


    emit updateHistogramValues(histogram);
    emit updateHistogramMaxValue(cumulation);
    emit updateTimeExpended(total_time);
    emit updateDistExpended(total_dist);
    emit updateNumSearchers( n_searchers );
    //emit updateNeeded();
    qApp->processEvents();
}

void Model::analyse()
{

    if (input_sets.size() != 3)
    {
        QErrorMessage errorMessage;
        errorMessage.showMessage("Supports only three datasets at a time.");
        errorMessage.exec();
        return;
    }

    bool generate_paths = true;

    if (space->getSearchType() == 0) generate_paths = false;

    //cout << "Thread " << pthread_self() << " creating new worker" << endl;

    workerThread = new QThread();
    worker = new Worker();
    worker->setGeneratePaths(generate_paths);
    worker->setModel(this);

    worker->moveToThread(workerThread);

         connect(workerThread, SIGNAL(started()), worker, SLOT(startWork()));
         connect(worker, SIGNAL(workFinished()), workerThread, SLOT(quit()), Qt::DirectConnection);
         connect(this, SIGNAL(wakeupWorker()), worker, SLOT(wakeup()));
         //connect(this, &Model::startWorker, worker, &Worker::doWork);
         //connect(worker, &Worker::workFinished(), this, &Model::handleWorkerFinish());
         workerThread->start();
         //emit startWorker();

         //cout << "Finished creating worker." << endl;
}

void Model::doAnalysis()
{
    // Calculate target interactions

    //cout << "Do Analysis" << endl;

    int count = 0;
    int id = 0; // keep track of the poistion
    int file_count = 0;

    vector< vector<double>* > distances_used;
    vector< vector<double>* > times_used;
    vector< vector<double>* > targets_found;
    vector< vector<double>* > time_efficiencies;
    vector< vector<double>* > dist_efficiencies;
    vector< vector<double>* > hopkins_values;
    QStringList dataset_labels;

    for (int i = 0; i < 15; i++)
    {
        distances_used.push_back(new vector<double>());
        times_used.push_back(new vector<double>());
        targets_found.push_back(new vector<double>());
        time_efficiencies.push_back(new vector<double>());
        dist_efficiencies.push_back(new vector<double>());
        hopkins_values.push_back(new vector<double>());
    }

    // h < 3
    for (int h = 0; h < input_sets.size(); h++)
    {
        int setid = h;

        QStringList* queued_input_files = input_sets[h];

        QStringList parts = queued_input_files->at(0).split(QDir::separator());
        dataset_labels.push_back(parts.at(parts.size()-2));


        for (int k = 0; k < queued_input_files->size(); k++)
        {
            int file_progress = 0;
            emit updateInputFileProcessing(file_count);

            cout << "Adding searchers from file " << queued_input_files->at(k).toStdString() << "..." << flush;
            addSearchers(queued_input_files->at(k));
            cout << "done" << endl;

            //cout << "Processing file " << file_count << endl;

            //for (int m = 0; m < 3; m++)
            {

                //  setSearchType(m);

                bool generate_paths = true;

                if (space->getSearchType() == 0) generate_paths = false;

                if (gui_ready)
                    if (generate_paths)
                    {
                        paths_mutex = true;
                        while ( !paths_mutex_seen )
                        {
                            //cout << "Checked path mutex" << endl;

                            //emit updateNeeded();
                            usleep(10);
                        }

                        clearPaths();
                        paths_mutex = false;
                        paths_mutex_seen = false;
                        //emit updateNeeded();

                        generate();

                        //cout << "clearPaths and generate new ones finished" << endl;
                    }

                // explorer different cluster radius values // 10, 20, 30, 40, 50
                for (int j = 0; j < 5; j++ )
                {

                    cluster_radius = j*10+10;
                    emit updateClusterRadius(cluster_radius);

                    int n_replicas = 100;

                    double dist_efficiency_mean;
                    double dist_efficiency_std;

                    for ( int i = 0; i < n_replicas; i++ )
                    {
                        qApp->processEvents();


                        distance_used = 0;
                        time_used = 0;
                        num_targets_found = 0;

                        placeTargets();

                        vector<Coordinate*> target_coords;

                        Target* tars = space->getTargets();

                        for (int ti = 0; ti < space->getNumTargets(); ti++)
                            target_coords.push_back(tars[ti].getCoordinates());

                        double hopkins_value = hopkins(target_coords);

                        int ntf = calcTargetsFound();

                        //sleep(1);

                        id = 5*setid+j+1;

                        int idx = id - 1;

                        distances_used[idx]->push_back(distance_used);
                        times_used[idx]->push_back(time_used);
                        targets_found[idx]->push_back(num_targets_found);
                        hopkins_values[idx]->push_back(hopkins_value);

                        count++;

                        time_efficiencies[idx]->push_back(num_targets_found/time_used);
                        dist_efficiencies[idx]->push_back(num_targets_found/distance_used);

                        double dist_mean = mean(*distances_used[idx]);
                        double time_mean = mean(*times_used[idx]);
                        double targets_found_mean = mean(*targets_found[idx]);
                        double dist_std = standardDeviation(*distances_used[idx]);
                        double time_std = standardDeviation(*times_used[idx]);
                        double targets_found_std = standardDeviation(*targets_found[idx]);
                        double time_efficiency_mean = mean(*time_efficiencies[idx]);
                        dist_efficiency_mean = mean(*dist_efficiencies[idx]);
                        double time_efficiency_std = standardDeviation(*time_efficiencies[idx]);
                        dist_efficiency_std = standardDeviation(*dist_efficiencies[idx]);

                        emit updateProgressBar((++file_progress)*100/(n_replicas*5));
                        //emit updateNeeded();
                        emit updateDistEfficiencyMean( dist_efficiency_mean ); //sum(targets_found)/sum(distances_used) );
                        emit updateTimeEfficiencyMean( time_efficiency_mean ); //sum(targets_found)/sum(times_used) );
                        emit updateDistEfficiencyStd( dist_efficiency_std );
                        emit updateTimeEfficiencyStd( time_efficiency_std );
                        qApp->processEvents();

                        FloatVector time_eff_sample;
                        FloatVector dist_eff_sample;

                        double time_eff_median = median(*time_efficiencies[idx]);
                        double time_eff_lower_quartile = lowerQuartile(*time_efficiencies[idx]);
                        double time_eff_upper_quartile = upperQuartile(*time_efficiencies[idx]);
                        double time_eff_min = min(*time_efficiencies[idx]);
                        double time_eff_max = max(*time_efficiencies[idx]);
                        double time_eff_interquartile_range = time_eff_upper_quartile - time_eff_lower_quartile;
                        double time_eff_inner_fence = 1.5*time_eff_interquartile_range;

                        vector<double> time_eff_outliers = calcOutliers(*time_efficiencies[idx], time_eff_inner_fence, time_eff_median);

                        double hopkins_mean = mean(*hopkins_values[idx]);

                        double dist_eff_median = median(*dist_efficiencies[idx]);
                        double dist_eff_lower_quartile = lowerQuartile(*dist_efficiencies[idx]);
                        double dist_eff_upper_quartile = upperQuartile(*dist_efficiencies[idx]);
                        double dist_eff_min = min(*dist_efficiencies[idx]);
                        double dist_eff_max = max(*dist_efficiencies[idx]);
                        double dist_eff_interquartile_range = dist_eff_upper_quartile - dist_eff_lower_quartile;
                        double dist_eff_inner_fence = 1.5*dist_eff_interquartile_range;

                        vector<double> dist_eff_outliers = calcOutliers(*dist_efficiencies[idx], dist_eff_inner_fence, dist_eff_median);

                        //cout << "{" << pthread_self() << "} Count: " << count << endl;
                        //cout << "Search Type: " << space->getSearchType() << endl;
                        // cout << "ID: " << id << endl;
                        //cout << "Replica: " << i << endl;
                        //cout << "Time Efficiencies vector: ";
                        //printVectorContents(time_efficiencies);
                        ////cout << "Cluster radius: " << j << endl;
                        ////cout << "Median" << time_eff_median << endl;
                        ////cout << "Upper Quartile: " << time_eff_upper_quartile << endl;
                        ////cout << "Lower Quartile: " << time_eff_lower_quartile << endl;
                        //cout << "Inner Fence: +-" << time_eff_inner_fence << endl;
                        //cout << "Max: " << time_eff_max << endl;
                        //cout << "Min: " << time_eff_min << endl;

                        if (time_efficiencies[idx]->size() > 4) // Make sure we have enough samples for which it makes sense to ask about quartiles
                        {
                            time_eff_sample.push_back(id);
                            time_eff_sample.push_back(setid);

                            time_eff_sample.push_back(time_eff_median);
                            time_eff_sample.push_back(time_eff_lower_quartile);
                            time_eff_sample.push_back(time_eff_upper_quartile);
                            //time_eff_sample.push_back(time_eff_min);
                            //time_eff_sample.push_back(time_eff_max);
                            time_eff_sample.push_back(time_eff_median-time_eff_inner_fence);
                            time_eff_sample.push_back(time_eff_median+time_eff_inner_fence);

                            //time_eff_sample.insert(time_eff_sample.end(), time_eff_outliers.begin(), time_eff_outliers.end());
                            time_eff_sample.push_back(hopkins_mean);
                            time_eff_sample.push_back(cluster_radius);

                            for (vector<double>::iterator it = time_eff_outliers.begin(); it != time_eff_outliers.end(); ++it)
                                time_eff_sample.push_back(*it);

                            //cout << "Outliers [" << time_eff_median-time_eff_inner_fence << "," << time_eff_median+time_eff_inner_fence << "]: "; printVectorContents(time_eff_outliers); cout << endl;

                            dist_eff_sample.push_back(id);
                            dist_eff_sample.push_back(setid);
                            dist_eff_sample.push_back(dist_eff_median);
                            dist_eff_sample.push_back(dist_eff_lower_quartile);
                            dist_eff_sample.push_back(dist_eff_upper_quartile);
                            //dist_eff_sample.push_back(dist_eff_min);
                            //dist_eff_sample.push_back(dist_eff_max);
                            dist_eff_sample.push_back(dist_eff_median-dist_eff_inner_fence);
                            dist_eff_sample.push_back(dist_eff_median+dist_eff_inner_fence);

                            dist_eff_sample.push_back(hopkins_mean);
                            dist_eff_sample.push_back(cluster_radius);

                            //dist_eff_sample.insert(dist_eff_sample.end(), dist_eff_outliers.begin(), dist_eff_outliers.end());

                            //time_eff_sample.insert(time_eff_sample.end(), time_eff_outliers.begin(), time_eff_outliers.end());
                            for (vector<double>::iterator it = dist_eff_outliers.begin(); it != dist_eff_outliers.end(); ++it)
                                dist_eff_sample.push_back(*it);

                            emit updateTimeEffSample(time_eff_sample);
                            emit updateDistEffSample(dist_eff_sample);
                            qApp->processEvents();
                        }

                    }

                    //cout << "ID: " << id << endl;

                    // for the result plot.
                    //cout << "Mean: " << dist_efficiency_mean << endl;
                    //cout << "STD: " << dist_efficiency_std << endl;
                }
                emit updateInputFileProcessed(file_count);
                qApp->processEvents();
                file_count++;

            }
        }
    }


    for (int i = 0; i < 5; i++)
    {
    cout << " distances used: "; printVectorContents(*distances_used[i]); cout << endl << endl;
    cout << " times used: "; printVectorContents(*times_used[i]); cout << endl << endl;
    cout << " targets found: "; printVectorContents(*targets_found[i]); cout << endl << endl;
    cout << " time_efficiencies: ";  printVectorContents(*time_efficiencies[i]); cout << endl << endl;
    cout << " distances efficiencies: ";  printVectorContents(*dist_efficiencies[i]); cout << endl << endl;
    cout << " hopkins values: ";  printVectorContents(*hopkins_values[i]); cout << endl << endl;
    }



    QDebugStream qout(std::cout, window->myTextEdit);

    double Sm1t = mean(*time_efficiencies[0]);
    double Sm2t = mean(*time_efficiencies[5]);
    double Sm3t = mean(*time_efficiencies[10]);

    double Sd1t = standardDeviation(*time_efficiencies[0]);
    double Sd2t = standardDeviation(*time_efficiencies[5]);
    double Sd3t = standardDeviation(*time_efficiencies[10]);

    double Sn1t = time_efficiencies[0]->size();
    double Sn2t = time_efficiencies[5]->size();
    double Sn3t = time_efficiencies[10]->size();

    double Sm1d = mean(*dist_efficiencies[0]);
    double Sm2d = mean(*dist_efficiencies[5]);
    double Sm3d = mean(*dist_efficiencies[10]);

    double Sd1d = standardDeviation(*dist_efficiencies[0]);
    double Sd2d = standardDeviation(*dist_efficiencies[5]);
    double Sd3d = standardDeviation(*dist_efficiencies[10]);


    double Sn1d = dist_efficiencies[0]->size();
    double Sn2d = dist_efficiencies[5]->size();
    double Sn3d = dist_efficiencies[10]->size();

    string S1_label = dataset_labels.at(0).toStdString();
    string S2_label = dataset_labels.at(1).toStdString();
    string S3_label = dataset_labels.at(2).toStdString();



    double alpha = 0.001;

    cout << "Student's T-Test Time Efficency " << S1_label << " and " << S2_label << endl;
    studentsTtestUsingBoost( S1_label, Sm1t, Sd1t, Sn1t, S2_label, Sm2t, Sd2t, Sn2t, alpha);
    studentsTtestDifferingSTDUsingBoost(S1_label, Sm1t, Sd1t, Sn1t, S2_label, Sm2t, Sd2t, Sn2t, alpha);

    cout << "Student's T-Test Distance Efficency " << S1_label << " and " << S2_label << endl;
    studentsTtestUsingBoost(S1_label, Sm1d, Sd1d, Sn1d, S2_label, Sm2d, Sd2d, Sn2d, alpha);
    studentsTtestDifferingSTDUsingBoost(S1_label, Sm1d, Sd1d, Sn1d, S2_label, Sm2d, Sd2d, Sn2d, alpha);

    cout << "Student's T-Test Time Efficency " << S1_label << " and " << S3_label << endl;
    studentsTtestUsingBoost(S1_label, Sm1t, Sd1t, Sn1t, S3_label, Sm3t, Sd3t, Sn3t, alpha);
    studentsTtestDifferingSTDUsingBoost(S1_label, Sm1t, Sd1t, Sn1t, S3_label, Sm3t, Sd3t, Sn3t, alpha);

    cout << "Student's T-Test Distance Efficency " << S1_label << " and " << S3_label << endl;
    studentsTtestUsingBoost(S1_label, Sm1d, Sd1d, Sn1d, S3_label, Sm3d, Sd3d, Sn3d, alpha);
    studentsTtestDifferingSTDUsingBoost(S1_label, Sm1d, Sd1d, Sn1d, S3_label, Sm3d, Sd3d, Sn3d, alpha);

    double p_twotail = 0;
    double p_lefttail = 0;
    double p_righttail = 0;

    real_1d_array sample1;
    real_1d_array sample2;
    real_1d_array sample3;


    alglib::ae_int_t sample1_size = (*time_efficiencies[0]).size();
    alglib::ae_int_t sample2_size = (*time_efficiencies[5]).size();
    alglib::ae_int_t sample3_size = (*time_efficiencies[10]).size();

    sample1.setcontent(sample1_size, (*time_efficiencies[0]).data());
    sample2.setcontent(sample2_size, (*time_efficiencies[5]).data());
    sample3.setcontent(sample3_size, (*time_efficiencies[10]).data());


    alglib::mannwhitneyutest(sample1, sample1_size, sample2, sample2_size, p_twotail, p_lefttail, p_righttail);

    cout <<
       "_______________________________________________\n"
       "Mann-Whitney test for two samples (does not assume populations are normally distributed)\n"
       "_______________________________________________\n\n";

    cout << S1_label << " vs. " << S2_label << endl;
    cout << "Twotail p-value: " << p_twotail << endl; //<< ". " << (p_twotail < alpha ? "REJECTED": "NOT-REJECTED") << endl;
    cout << "Righttail p-value: " << p_righttail  << endl; //<< (p_righttail < alpha ? "REJECTED": "NOT-REJECTED") << endl;
    cout << "Lefttail p-value: " << p_lefttail  << endl; //<< (p_lefttail < alpha ? "REJECTED": "NOT-REJECTED") << endl;

    alglib::mannwhitneyutest(sample1, sample1_size, sample3, sample3_size, p_twotail, p_lefttail, p_righttail);

    cout << "_______________________________________________\n" << endl;
    cout << S1_label << " vs. " << S3_label << endl;
    cout << "Twotail p-value: " << p_twotail  << endl; //<< ". " << (p_twotail < alpha ? "REJECTED": "NOT-REJECTED") << endl;
    cout << "Righttail p-value: " << p_righttail  << endl; //<< (p_righttail < alpha ? "REJECTED": "NOT-REJECTED") << endl;
    cout << "Lefttail p-value: " << p_lefttail  << endl; //<< (p_lefttail < alpha ? "REJECTED": "NOT-REJECTED") << endl;
    cout << "_______________________________________________\n" << endl;

}

int Model::calcTargetsFound()
{

    Target* targets = space->getTargets();
    // Check whether searchers have found targets
      for ( int j = 0; j < space->getNumTargets(); j++ )
      {
        targets[j].setNotFound();
      }


    // Calculate target interactions
    int n_searchers = space->getNumSearchers();
    Searcher* searchers = space->getSearchers();

    for ( int i = 0; i < n_searchers; i++ )
    {
        Searcher s = searchers[i];
        vector<Coordinate*> path = s.getPath();

        if (path.size() >= 2)
        for (int j = 0; j < path.size()-1; j++)
        {

   //qApp->processEvents();
                ////cout << "j = " << j << " of " << path.size() << endl;
                Coordinate* start = path[j];
                Coordinate* end = path[j+1];
                space->testpath( start->getX(), start->getY(), start->getZ(), end->getX(), end->getY(), end->getZ() );
        }

    }

    int n_targets_found = 0;
    // Check whether searchers have found targets
      for ( int j = 0; j < space->getNumTargets(); j++ )
      {
          if (targets_with_replacement)
          n_targets_found += targets[j].getTotalContacts();
            else
          if (targets[j].isFound()) n_targets_found++;
      }

    double d = this->getSearchSpace()->getTotalDistanceTravelled();
    double t = this->getSearchSpace()->getTotalTimeTravelled();

    t *= time_resolution;

    d == 0 ? dist_efficiency = 0: dist_efficiency = n_targets_found/d;
    t == 0 ? time_efficiency = 0: time_efficiency = n_targets_found/t;
    distance_used = d;
    time_used = t;
    num_targets_found = n_targets_found;

    /*
    //cout << "Time travelled: " << t << endl;
    //cout << "Distance travelled: " << d << endl;

    //cout << "Time Efficiency: " << time_efficiency << endl;
    //cout << "Distance Efficieny: " << dist_efficiency << endl;
*/

    //emit updateNeeded();
    return n_targets_found;
}

 float Model::getDetectionRadius()
 {
     return detection_radius;
 }

 int Model::getNumClusters()
 {
     return n_clusters;
 }

 int Model::getNumTargetsPerCluster()
 {
     return targets_per_cluster;
 }

 int Model::getNumSearchers()
 {
     return this->getSearchSpace()->getNumSearchers();
 }

 int Model::getNumTargets()
 {
     return this->getSearchSpace()->getNumTargets();
 }

 void Model::generateSearch()
 {
     //cout << "Model:generateSearch() called" << endl;

       }

 void Model::clearPaths()
 {
     pthread_t id = pthread_self();
     //cout << "Model:clearPaths() called by thread " << id << endl;
     int n_searchers = space->getNumSearchers();
     Searcher* searchers = space->getSearchers();
     for (int i = 0; i < n_searchers; i++) searchers[i].clearPath();
 }

 void Model::repeat()
 {
    removeTargets();
    placeTargets();
    emit updateNeeded();
 }

 void Model::placeTargets()
 {    // Add targets
     // Create targets

     //cout << "Model: Targets being placed" << endl;

     float x, y, z;

     //cout << "Num Targets: " << max_targets << endl;
     //cout << "Num Clusters: " << n_clusters << endl;

     //removeTargets();

     int n_targets = 200*scale_number_of_targets;
     //Target* targets = getSearchSpace()->getTargets();
     Target* targets = new Target[n_targets];
     getSearchSpace()->setTargets( targets, n_targets );

     emit updateDensity(n_targets/getSearchSpace()->getVolume());

   //  getSearchSpace()->setNumTargets(n_targets);

     n_clusters = n_targets/targets_per_cluster+1; // +1 for overflow cluster
     emit updateNumClusters(n_clusters);

    // cout << "--------------------------" << endl;
    // cout << "n_targets: " << n_targets << endl;
    // cout << "targets_per_cluster: " << targets_per_cluster << endl;
    // cout << "n_clusters: " << n_clusters << endl;
    // cout << "--------------------------" << endl;


 // Initialize targets
     for ( int j = 0; j < space->getNumTargets(); j++ )
       {
         targets[j].setNotFound();
       }

     int x_min_bound = ceil(getSearchSpace()->getXMinBound());
     int y_min_bound = ceil(getSearchSpace()->getYMinBound());
     int z_min_bound = ceil(getSearchSpace()->getZMinBound());

     int x_max_bound = ceil(getSearchSpace()->getXMaxBound());
     int y_max_bound = ceil(getSearchSpace()->getYMaxBound());
     int z_max_bound = ceil(getSearchSpace()->getZMaxBound());

     int center_x = 0;
     int center_y = 0;
     int center_z = 0;

     int n_placed = 0;

     for (int i = 0; i < n_clusters; i++)
     {                     
             center_x = rand()%(x_min_bound-x_max_bound)+x_min_bound;
             center_y = rand()%(y_min_bound-y_max_bound)+y_min_bound;
             center_z = rand()%(z_min_bound-z_max_bound)+z_min_bound;

         //cout << "Cluster Center: " << center_x << ", " <<  center_y << ", " <<  center_z << endl;

         int j = 0;
         while( j < targets_per_cluster)
         {
                 x = (rand()%2000*cluster_radius/1000.0-cluster_radius);
                 y = (rand()%2000*cluster_radius/1000.0-cluster_radius);
                 z = (rand()%2000*cluster_radius/1000.0-cluster_radius);


              // cout << "Potential: " << x << ", " << y << ", " << z << endl;


             if ( sqrt(x*x + y*y + z*z) < cluster_radius )
             {
                 if (n_placed >= n_targets) break;

                 x += center_x;
                 y += center_y;
                 z += center_z;

                 if ( x < x_max_bound && x > x_min_bound
                 && y < y_max_bound && y > y_min_bound
                 && z < z_max_bound && z > z_min_bound )
                 {
                 //cout << "scaled " << x << ", " << y << ", " << z << endl;

                 targets[n_placed].setXPos(x);
                 targets[n_placed].setYPos(y);
                 targets[n_placed].setZPos(z);
                 j++;

                 n_placed++;
                 //cout << "Placed target in cluster" << "(" << n_placed << ") " << " at "<< "(" << x << ", " << y << ", " << z << ")" <<endl;
                 }
             }
         }
     }


     // cout << n_placed << ", " << n_targets << endl;

     emit updateNumTargets(n_placed);

     /*
     cout << "Successfully placed " << n_placed << " targets" << endl;


     cout << "x min bound: " << x_min_bound << endl;
     cout << "x max bound: " << x_max_bound << endl;

     cout << "y min bound: " << y_min_bound << endl;
     cout << "y max bound: " << y_max_bound << endl;

     cout << "z min bound: " << z_min_bound << endl;
     cout << "z max bound: " << z_max_bound << endl;
     */


 }

 void Model::removeTargets()
 {
     targets_mutex = true;
    // emit updateNeeded();

     if ( gui_ready )
     {
     while ( !targets_mutex_seen )
     {
         //cout << "Checked target mutex" << endl;
         //emit updateNeeded();
         qApp->processEvents();

     }
     }

     if (space->getTargets() != NULL)
     {
        space->clearTargets();
        //delete [] space->getTargets();
     }


     targets_mutex = false;
     targets_mutex_seen = false;
 }

 void Model::removeSearchers()
 {
     searchers_mutex = true;
     //emit updateNeeded();

     if ( gui_ready )
     {
     while ( !searchers_mutex_seen )
     {
         //cout << "Checked searcher mutex" << endl;
         //emit updateNeeded();
         qApp->processEvents();
         usleep(10);
     }
     }

     if (space->getSearchers() != NULL)
     {
        space->clearSearchers();
        //delete [] space->getTargets();
     }


     searchers_mutex = false;
     searchers_mutex_seen = false;
 }

 void Model::handleWorkerFinish()
 {
     cout << "Worker finish being handled by thread " << pthread_self() << endl;

     if (workerThread->isRunning())
     {
         cout << "Worker thread " << " is still running." << endl;
         workerThread->exit();
        // workerThread.wait();
     }

     // cout << "Worker thread finished. Deleting." << endl;
     //workerThread.deleteLater();

     calcTargetsFound();

 }

 void Model::placeSearchers()
 {    // Add targets
     // Create targets

     cout << "Model::placeSearchers() called." << endl;

     float x, y, z;
     int n_searchers = 50;
     Searcher* searchers = new Searcher[n_searchers];

     int x_min_bound = ceil(getSearchSpace()->getXMinBound());
     int y_min_bound = ceil(getSearchSpace()->getYMinBound());
     int z_min_bound = ceil(getSearchSpace()->getZMinBound());

     int x_max_bound = ceil(getSearchSpace()->getXMaxBound());
     int y_max_bound = ceil(getSearchSpace()->getYMaxBound());
     int z_max_bound = ceil(getSearchSpace()->getZMaxBound());

     // create searchers
     for (int i = 0; i < n_searchers; i++)
     {
         x = rand()%(x_min_bound-x_max_bound)+x_min_bound;
         y = rand()%(y_min_bound-y_max_bound)+y_min_bound;
         z = rand()%(z_min_bound-z_max_bound)+z_min_bound;

         searchers[i].setXPos(x);
         searchers[i].setYPos(y);
         searchers[i].setZPos(z);
     }

     getSearchSpace()->setSearchers(searchers, n_searchers);

     //emit updateNeeded();
 }

void Model::setSearchType( int v )
{
    space->setSearchType( v );
}

void Model::saveSettings()
{
    // Write settings to powersearch.cfg

    string path = app_path + "powersearch.cfg";
    cout << "Model: saving settings to " << path << endl;

    ofstream config_file;
    config_file.open(path.c_str());

    config_file << working_path << endl;

    cout << "Saved working path value " << working_path << endl;

    config_file.close();
}

void Model::readSettings()
{
    // Read settings from powersearch.cfg

    string path = app_path + "powersearch.cfg";

    cout << "Model: Reading settings from " << path << endl;

    ifstream config_file;
    config_file.open(path.c_str());

    getline(config_file, working_path);
    config_file.close();

    cout << "Model: working_path set to " << working_path << endl;
}

string Model::getWorkingPath()
{
    return working_path;
}

double Model::sum( vector<double> values )
{
        double vsum = std::accumulate( values.begin(), values.end(), 0.0f );
        return vsum;
}

void Model::printVectorContents( vector<double> values )
{
    for(std::vector<double>::iterator j=values.begin(); j!=values.end(); ++j)
        cout << *j << " ";
    cout << endl;
}

double Model::mean( vector<double> values )
{
    if (values.size() == 0) return 0;

    double vsum = sum(values);
    return vsum/values.size();
}

double Model::median(vector<double> values)
{
    if (values.size() == 0) return -99999999;
    if (values.size() == 1) return values[0];

    double m = 0;
    sort(values.begin(), values.end());

    if (values.size()%2 != 0) m  = values[values.size()/2]; // odd case
    else m = 0.5*(values[values.size()/2]+values[values.size()/2-1]); // even case

    return m;
}

double Model::lowerQuartile(vector<double> values)
{

    sort(values.begin(), values.end());
    int median_index = -1;
    if (values.size()%2 != 0) median_index  = values.size()/2;
    else median_index = values.size()/2-1;
    vector<double> nv(values.begin(), values.begin()+median_index);

    double m = median(nv);

    return m;
}

double Model::upperQuartile(vector<double> values)
{
    sort(values.begin(), values.end());
    int median_index = -1;
    if (values.size()%2 != 0) median_index = values.size()/2;
    else median_index = values.size()/2-1;

    vector<double> nv(values.begin()+median_index, values.end() );

    double m = median(nv);

    return m;

}

double Model::min(vector<double> values)
{
    double min = 100000;

    for(std::vector<double>::iterator j=values.begin(); j!=values.end(); ++j) if (*j < min) min = *j;

    return min;
}

double Model::max(vector<double> values)
{

    double max = -100000;

    for(std::vector<double>::iterator j=values.begin(); j!=values.end(); ++j) if (*j > max) max = *j;

    return max;
}

double Model::variance( vector<double> values )
{
    if (values.size() == 0) return 0;

    double m = mean(values);

    vector<double> terms;

    for(std::vector<double>::iterator it = values.begin(); it != values.end(); ++it)
    {
        double term = (*it - m);
        terms.push_back(term*term);
    }

    return sum(terms)/terms.size();
}

double Model::standardDeviation( vector<double> values )
{
    return sqrt(variance(values));
}

void Model::paintingGLFinished()
{
    cout << "Model: paintingGLFinished called" << endl;
    emit wakeupWorker();
}

void Model::queueInputFiles(QStringList l)
{
    //emit clearInputFilesToProcess();

    QString fp = *(l.begin());

    std::string filepath = fp.toUtf8().constData();

    //cout << "Reading file " << filepath << endl;

    QFileInfo info(QString::fromStdString(filepath));

    working_path = info.dir().absolutePath().toStdString();

    cout << "Model: working_path set to " << working_path << endl;


    // Only allow three input sets for now
    if ( input_sets.size() > 2)
    {
        QErrorMessage errorMessage;
        errorMessage.showMessage("Supports only three datasets at a time.");
        errorMessage.exec();
        return;
    }

    emit addInputFilesToProcess(l);


    QStringList* lptr = new QStringList(l);
    input_sets.push_back(lptr);


}


void Model::setTargetsWithReplacement(bool v)
{
    targets_with_replacement = v;
    cout << "With replacement set to ";
    if (targets_with_replacement) cout << "\"yes\"" << endl;
    else cout << "\"no\"" << endl;
}

void Model::addInputSet(QStringList* input_file_list)
{
    input_sets.push_back(input_file_list);
}

vector<double> Model::calcOutliers(vector<double> values, double inner_fence, double median)
{
    vector<double> outliers;

    for(std::vector<double>::iterator j=values.begin(); j!=values.end(); ++j)
        if (*j < median-inner_fence || *j > median+inner_fence) outliers.push_back(*j);

    return outliers;
}

double Model::hopkins(vector<Coordinate*> S)
{
    // The ideas is to compare the sum of the nearest neighbor distances for the locations of interest with the nearest neighbor distances of a uniformly distributed set.

    // 1) Generate a vector of uniformly distributed points in the search space, R.
    // 2) For each coordinate in S find the nearest neighbor in R.
    // 3) Keep a running sum of these nearest neighbor distances, W.
    // 4) For each coordinate in S find the nearest neighbor in locations.
    // 5) keep a running sum of these nearest neighbor distances, U.
    // 6) The normalized Hopkins statistic is U/(U+W)
    // 7) Values close to 1/2 are not clustered, values close to 0 are more clustered.

    vector<Coordinate*> R;

    double U = 0;
    double W = 0; // Running sum for uniformly distributed nearst neighbor distances

    double min_x = space->getXMinBound();
    double max_x = space->getXMaxBound();

    double min_y = space->getYMinBound();
    double max_y = space->getYMaxBound();

    double min_z = space->getZMinBound();
    double max_z = space->getZMaxBound();

    for(int i = 0; i < S.size(); i++)
    {
        // Generate random coordinates
        double x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
        double y = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
        double z = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);

        // Scale x to be between the limits of the field
        x = x*(max_x-min_x)+min_x;
        y = y*(max_y-min_y)+min_y;
        z = z*(max_z-min_z)+min_z;

        double time = 0;

        Coordinate* c = new Coordinate(x, y, z, time);
        R.push_back(c);
    }

    // Compare nearest neighbor distance sums for the uniform and observed distributions
    for (std::vector<Coordinate*>::iterator Sit = S.begin() ; Sit != S.end(); ++Sit)
    {
        double min_R_distance = 10000000;
        double min_S_distance = 10000000;

        for (std::vector<Coordinate*>::iterator Rit = R.begin() ; Rit != R.end(); ++Rit)
        {
            double dist = space->norm((*Sit)->getX()-(*Rit)->getX(),(*Sit)->getY()-(*Rit)->getY(),(*Sit)->getZ()-(*Rit)->getZ());
            if (dist < min_R_distance)
            {
                min_R_distance = dist;
            }
        }

        for (std::vector<Coordinate*>::iterator Sit2 = S.begin() ; Sit2 != S.end(); ++Sit2)
        {
            double dist = space->norm((*Sit)->getX()-(*Sit2)->getX(),(*Sit)->getY()-(*Sit2)->getY(),(*Sit)->getZ()-(*Sit2)->getZ());
            if ((*Sit) != (*Sit2) && dist < min_S_distance)
            {
                min_S_distance = dist;
            }
        }

        // Add distance to the nearest neighbor to the sum
        W += min_R_distance;
        U += min_S_distance;
    }

    double result = U/(U+W);

    //cout << "Hopkins Clustering Statistic: U//(U+W) = " << U << "//" << "(" << U << "+" << W << ") = " << result << endl;

    return result;

}

vector<Coordinate*> Model::getHopkinsPoints()
{
    return hopkins_points;
}

// Assumes the distributions (in the null hypothesis?) have the same variance, does not assume the number of samples are equal
double Model::studentsTtest(vector<double> sample1, vector<double> sample2)
{
    // t is result of the t test.
    // S1S2 is the grand (or pooled) standard deviation
    // X_bars are the means of sample 1 and 2
    // S_x is the standard deviation for sample x
    // n_x is the number of samples in sample 1 and 2
    // Note: for significance testing the degrees of freedom are 2n-2.

    double n_1 = sample1.size();
    double n_2 = sample2.size();

    double X_bar_1 = mean(sample1);
    double X_bar_2 = mean(sample2);

    double S_1 = standardDeviation(sample1);
    double S_2 = standardDeviation(sample2);
    double S_1_squared = S_1*S_1;
    double S_2_squared = S_2*S_2;


    double S1S2 = sqrt(((n_1-1)*S_1_squared+(n_2-1)*S_2_squared)/(n_1+n_2-2));

    double t = (X_bar_1 = X_bar_2)/(S1S2*sqrt(1/n_1+1/n_2));

    return t;
}

void Model::studentsTtestUsingBoost(
        string S1_label,
        double Sm1,
        double Sd1,
        unsigned Sn1,
        string S2_label,
        double Sm2,
        double Sd2,
        unsigned Sn2,
        double alpha)
{
   //
   // Sm1 = Sample Mean 1.
   // Sd1 = Sample Standard Deviation 1.
   // Sn1 = Sample Size 1.
   // Sm2 = Sample Mean 2.
   // Sd2 = Sample Standard Deviation 2.
   // Sn2 = Sample Size 2.
   // alpha = Significance Level.
   //
   // A Students t test applied to two sets of data.
   // We are testing the null hypothesis that the two
   // samples have the same mean and that any difference
   // if due to chance.
   // See http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm
   //

   // Print header:
   cout <<
      "_______________________________________________\n"
      "Student t test for two samples (equal variances)\n"
      "_______________________________________________\n\n";
   cout << setprecision(5);
   cout << setw(55) << left << "Number of Observations ("+S1_label+")" << "=  " << Sn1 << "\n";
   cout << setw(55) << left << S1_label+" Mean" << "=  " << Sm1 << "\n";
   cout << setw(55) << left << S1_label+" Standard Deviation" << "=  " << Sd1 << "\n";
   cout << setw(55) << left << "Number of Observations ("+S2_label+")" << "=  " << Sn2 << "\n";
   cout << setw(55) << left << S2_label +" Mean" << "=  " << Sm2 << "\n";
   cout << setw(55) << left << S2_label +" Standard Deviation" << "=  " << Sd2 << "\n";
   //
   // Now we can calculate and output some stats:
   //
   // Degrees of freedom:
   double v = Sn1 + Sn2 - 2;
   cout << setw(55) << left << "Degrees of Freedom" << "=  " << v << "\n";
   // Pooled variance:
   double sp = sqrt(((Sn1-1) * Sd1 * Sd1 + (Sn2-1) * Sd2 * Sd2) / v);
   cout << setw(55) << left << "Pooled Standard Deviation" << "=  " << v << "\n";
   // t-statistic:
   double t_stat = (Sm1 - Sm2) / (sp * sqrt(1.0 / Sn1 + 1.0 / Sn2));
   cout << setw(55) << left << "T Statistic" << "=  " << t_stat << "\n";
   //
   // Define our distribution, and get the probability:
   //
   students_t dist(v);
   double q = cdf(complement(dist, fabs(t_stat)));
   cout << setw(55) << left << "Probability that difference is due to chance" << "=  "
      << setprecision(5) << scientific << 2 * q << "\n\n";
   //
   // Finally print out results of alternative hypothesis:
   //
   cout << setw(55) << left <<
      "Results for Alternative Hypothesis and alpha" << "=  "
      << setprecision(4) << fixed << alpha << "\n\n";
   cout << "Alternative Hypothesis          Conclusion\n";
   cout << S1_label << " Mean != "<<S2_label<<" Mean       " ;
   if(q < alpha / 2)
      cout << "NOT REJECTED\n";
   else
      cout << "REJECTED\n";
   cout << S1_label << " Mean <  "<<S2_label<<" Mean       ";
   if(cdf(dist, t_stat) < alpha)
      cout << "NOT REJECTED\n";
   else
      cout << "REJECTED\n";
   cout << S1_label << " Mean >  "<<S2_label<<" Mean       ";
   if(cdf(complement(dist, t_stat)) < alpha)
      cout << "NOT REJECTED\n";
   else
      cout << "REJECTED\n";
   cout << endl << endl;
}

void Model::studentsTtestDifferingSTDUsingBoost(
        string S1_label,
        double Sm1,
        double Sd1,
        unsigned Sn1,
        string S2_label,
        double Sm2,
        double Sd2,
        unsigned Sn2,
        double alpha)
{
   //
   // Sm1 = Sample Mean 1.
   // Sd1 = Sample Standard Deviation 1.
   // Sn1 = Sample Size 1.
   // Sm2 = Sample Mean 2.
   // Sd2 = Sample Standard Deviation 2.
   // Sn2 = Sample Size 2.
   // alpha = Significance Level.
   //
   // A Students t test applied to two sets of data.
   // We are testing the null hypothesis that the two
   // samples have the same mean and that any difference
   // if due to chance.
   // See http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm
   //

   // Print header:
   cout <<
      "_________________________________________________\n"
      "Student t test for two samples (unequal variances)\n"
      "_________________________________________________\n\n";
   cout << setprecision(5);



   cout << setw(55) << left << "Number of Observations ("+S1_label+")" << "=  " << Sn1 << "\n";
   cout << setw(55) << left << S1_label + " Mean" << "=  " << Sm1 << "\n";
   cout << setw(55) << left << S1_label + " Standard Deviation" << "=  " << Sd1 << "\n";
   cout << setw(55) << left << "Number of Observations ("+S2_label+")" << "=  " << Sn2 << "\n";
   cout << setw(55) << left << S2_label+" Mean" << "=  " << Sm2 << "\n";
   cout << setw(55) << left << S2_label+" Standard Deviation" << "=  " << Sd2 << "\n";
   //
   // Now we can calculate and output some stats:
   //
   // Degrees of freedom:
   double v = Sd1 * Sd1 / Sn1 + Sd2 * Sd2 / Sn2;
   v *= v;
   double t1 = Sd1 * Sd1 / Sn1;
   t1 *= t1;
   t1 /=  (Sn1 - 1);
   double t2 = Sd2 * Sd2 / Sn2;
   t2 *= t2;
   t2 /= (Sn2 - 1);
   v /= (t1 + t2);
   cout << setw(55) << left << "Degrees of Freedom" << "=  " << v << "\n";
   // t-statistic:
   double t_stat = (Sm1 - Sm2) / sqrt(Sd1 * Sd1 / Sn1 + Sd2 * Sd2 / Sn2);
   cout << setw(55) << left << "T Statistic" << "=  " << t_stat << "\n";
   //
   // Define our distribution, and get the probability:
   //
   students_t dist(v);
   double q = cdf(complement(dist, fabs(t_stat)));
   cout << setw(55) << left << "Probability that difference is due to chance" << "=  "
      << setprecision(3) << scientific << 2 * q << "\n\n";
   //
   // Finally print out results of alternative hypothesis:
   //
   cout << setw(55) << left <<
      "Results for Alternative Hypothesis and alpha" << "=  "
      << setprecision(4) << fixed << alpha << "\n\n";
   cout << "Alternative Hypothesis              Conclusion\n";
   cout << S1_label << " Mean != "<<S2_label<<" Mean       " ;
   if(q < alpha / 2)
      cout << "NOT REJECTED\n";
   else
      cout << "REJECTED\n";
   cout << S1_label << " Mean <  "<<S2_label<<" Mean       ";
   if(cdf(dist, t_stat) < alpha)
      cout << "NOT REJECTED\n";
   else
      cout << "REJECTED\n";
   cout << S1_label << " Mean >  "<<S2_label<<" Mean       ";
   if(cdf(complement(dist, t_stat)) < alpha)
      cout << "NOT REJECTED\n";
   else
      cout << "REJECTED\n";
   cout << endl << endl;
}

