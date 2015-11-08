#include "Searcher.h"
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <boost/math/special_functions/fpclassify.hpp>

using namespace std;

Searcher::Searcher()
{
    //cout << "Searcher: Searcher created" << endl;
    first_contact_time = 1000000.0f;
}

void Searcher::addToPath(Coordinate* c)
{
    //cout << "Position changed: id=" << getID() << " <" << getXPos() << ", " << getYPos() << ", " << getZPos() << ">" << endl;

   if (c == NULL)
   {
       cout << "Tried to add a NULL coordinate to the path. Skipping." << endl;
       return;
   }

    path.push_back(c);

}

void Searcher::scalePath(float min_x, float max_x, float min_y, float max_y, float min_z, float max_z)
{
    for (vector<Coordinate*>::iterator it = path.begin() ; it != path.end(); ++it)
    {
        (*it)->setX(2*(*it)->getX()/max_x-1);
        (*it)->setY(2*(*it)->getY()/max_y-1);
        (*it)->setZ(2*(*it)->getZ()/max_z-1);

        //cout << (*it)->getX() << " " << (*it)->getY() << " " << (*it)->getZ() << endl;
    }

 }

vector<Coordinate*> Searcher::getPath()
{
    return path;
}

void Searcher::setPath( vector<Coordinate*> v )
{
    path = v;
}

float Searcher::getDistanceTravelled()
{
    float dist = 0;
    if (path.size() > 2)
    for(vector<Coordinate*>::iterator it = path.begin(); it != path.end()-1; ++it)
    {
         dist += (*it)->getDistance(*(it+1));
        // if (boost::math::isnan(dist))
        // {
        //     cout << "Searcher::getDistanceTravelled(): dist is NaN" << endl;
        //     cout << "Distance between: <"<< (*it)->getX() <<","<< (*it)->getY() <<","<< (*it)->getZ() <<"> and <"<< (*it+1)->getX() <<","<< (*it+1)->getY() <<","<< (*it+1)->getZ() <<">" << endl;
        // }


    }
    else return 0;

    return dist;
}

float Searcher::getTimeTravelled()
{
    return path.size();
}

void Searcher::clearPath()
{
   //cout << "Searcher at address " << this << " is clearing path" << endl;
   path.clear();

   returnToStartingPosition();
   first_contact_time = 1000000.0f;

   //cout << "x:" << x_pos << " ";
   //cout << "y:" << y_pos << " ";
   //cout << "z:" << z_pos;
   //cout << endl;
}

void Searcher::addFoundTarget(Target* t, float time)
{
    if (first_contact_time > 1000000.0f-100)
    {
        first_contact_time = time;
    }

    if (!isTargetFound(t))
    {
        //cout << "Searcher: target added to searcher." << endl;
        found_targets.push_back(t);
        //cout << found_targets.size() << endl;
    }
}

bool Searcher::isTargetFound(Target* t)
{
    for(int i = 0; i < found_targets.size(); i++)
    {
        //cout << "found_targets[i]: " << found_targets[i] << endl;
        //cout << "t: " << t << endl;
        if (found_targets[i] == t) return true;
    }
    return false;
}

void Searcher::clearFoundTargets()
{
    first_contact_time = 1000000.0f;
    found_targets.clear();
}

int Searcher::getNumTargetsFound()
{
    //cout << "Searcher::getNumTargetsFound(): " << found_targets.size() << endl;
    return found_targets.size();
}

double Searcher::getFirstTargetContactTime()
{
    return first_contact_time;
}

float Searcher::getDistEfficiency()
{
    float dist_used = getDistanceTravelled();
    return found_targets.size()/dist_used;
}

float Searcher::getTimeEfficiency()
{
    float time_used = path.size();
    return found_targets.size()/time_used;
}

vector<double> Searcher::calcSpeeds()
{
    vector<Coordinate*> path = getPath();
    int path_length = path.size();

    return speeds;
}

vector<double> Searcher::calcAngles()
{

    vector<Coordinate*> path = getPath();
    int path_length = path.size();



    return angles;
}


vector<double> Searcher::getSpeeds()
{
    return speeds;
}

vector<double> Searcher::getAngles()
{
    return angles;
}
