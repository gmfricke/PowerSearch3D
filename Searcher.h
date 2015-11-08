#ifndef SEARCHER_H
#define SEARCHER_H

#include "Agent.h"
#include "Coordinate.h"
#include <vector>
#include <string>


using namespace std;

class Target;

class Searcher : public Agent
{
public:
    Searcher();
    Searcher(string path);
    void addToPath(Coordinate* c);
    vector<Coordinate*> getPath();
    void setPath( vector<Coordinate*> v );
    void scalePath(float min_x, float max_x, float min_y, float max_y, float min_z, float max_z);
    float getDistanceTravelled();
    float getTimeTravelled();
    void addFoundTarget(Target* t, float time);
    void clearPath();
    void clearFoundTargets();
    int getNumTargetsFound();
    bool isTargetFound(Target* t);
    double getFirstTargetContactTime();
    float getTimeEfficiency();
    float getDistEfficiency();
    vector<double> getSpeeds();
    vector<double> getAngles();

    vector<double> calcAngles();
    vector<double> calcSpeeds();

private:
    vector<Target*> found_targets;
    double first_contact_time;
    vector<double> angles;
    vector<double> speeds;

};

#endif // SEARCHER_H
