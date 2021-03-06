#ifndef SEARCHSPACE_H
#define SEARCHSPACE_H
#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include "Coordinate.h"

class Target;
class Searcher;
class Agent;

using namespace boost;
using namespace std;

class SearchSpace
{
    public:
    SearchSpace(int nt);
    Target* getTargets();
    Target* getTargetsCopy(int& n_targets);
    Searcher* getSearchersCopy(int& n_searchers);
    void setTargets(Target* targets, int n);
    Searcher* getSearchers();
    void setSearchers( Searcher* s, int n );
    vector<double> getSpeeds();
    vector<double> getAngles();
    void setBounds(float, float, float, float, float, float);

    float getBoundingHeight();
    float getBoundingDepth();
    float getBoundingWidth();

    float getXMinBound();
    float getYMinBound();
    float getZMinBound();

    float getXMaxBound();
    float getYMaxBound();
    float getZMaxBound();

    int getNumSearchers();
    int getNumTargets();
    void clearTargets();
    void clearSearchers();

    void setNumTargets(int v);

    double getTotalDistanceTravelled();
    double getTotalTimeTravelled();
    vector<Target*> testpath(float start_x, float start_y, float start_z, float end_x, float end_y, float end_z);
    float dot_prod(float x0, float y0, float z0, float x1, float y1, float z1);
    void cross_prod(float ux, float uy, float uz, float vx, float vy, float vz, float& rx, float& ry, float& rz);
    float norm(float x, float y, float z);
    float angle(float a_x, float a_y, float a_z, float b_x, float b_y, float b_z);
    void setTargetDetectRadius(float v);

    float getSpeedFromObserved();
    float getAngleFromObserved();
    float getFromGammaPDF();
    float getFromGammaPDF(float shape, float scale);
    float getFromNormalPDF();
    float getFromMaxwellPDF();
    void printMatrix(float** M);
    float getEuclideanDistance( Agent* a1, Agent* a2 );
    float getFromLogNormalPDF();
    float getFromNormalPDF(float mu, float sigma);
    float getFromGenParetoPDF();
    float getFromBetaPDF();
    float getFromPowerLawPDF();//double x_min, double x_max, double alpha);
    float getFromHarrisPowerLawPDF();//double x_min, double x_max, double alpha);
    float getFromExponentialPDF();
    float getFromObserved();
    float** matmul4by4(float A[4][4], float B[4][4]);
    Coordinate rotateAboutVector(Coordinate v, Coordinate axis, float angle);
    void setTimeResolution(float v);
    float getTimeResolution();
    void setCurrentTimeStep(int t);
    int getCurrentTimeStep();
    float getCurrentTime();

    void setObservedAngles(vector<double>);
    vector<double> getObservedAngles();
    void setObservedSpeeds(vector<double>);
    std::vector<double> getObservedSpeeds();

    void setSearchType(int v);
    int getSearchType();
    string getSearchTypeName();

    vector<double> getFirstContactTimes();

    virtual float getVolume(){}

    ~SearchSpace();

protected:
    Target* targets;
    Searcher* searchers;
    int n_searchers;
    int n_targets;
    double total_distance_travelled;
    float detection_radius;
    int search_pattern;  
    float x_min_bound;
    float y_min_bound;
    float z_min_bound;
    float x_max_bound;
    float y_max_bound;
    float z_max_bound;
    float current_time;
    int search_type;
    float time_resolution;
    float current_time_step;
    vector<double> observed_speeds;
    vector<double> observed_angles;

private:
    mt19937 rng;

};

#endif // SEARCHSPACE_H
