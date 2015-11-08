#include "Agent.h"
#include "Coordinate.h"
#include <cmath>
#include <stdlib.h>
#include <iostream>

using namespace std;

static int agent_count = 0;

Agent::Agent()
{
    this->id = agent_count++;
    speed = 0.1;
        
    starting_time = 0;
    starting_x = 0;
    starting_y = 0;
    starting_z = 0;

	radius = 0.01;

    red = rand()%128+127;
    blue = rand()%128+127;
    green = rand()%128+127;

    //red = 255;
    //blue = 255;
    //green = 255;
}

int Agent::getID()
{
    return id;
}

float Agent::getRadius()
{
    return radius;
}

void Agent::setRadius(float r)
{
    radius = r;
}

void Agent::setXPos( float x )
{
    x_pos = x;
}

void Agent::setYPos( float y )
{
    y_pos = y;
}

void Agent::setZPos( float z )
{
    z_pos = z;
}

void Agent::setTime( float t )
{
    time = t;
}

void Agent::setStartingXPos( float x )
{
    starting_x = x;
}

void Agent::setStartingYPos( float y )
{
    starting_y = y;
}

void Agent::setStartingZPos( float z )
{
    starting_z = z;
}

void Agent::setStartingTime( float t )
{
    starting_time = t;
}

float Agent::getXPos()
{
    return x_pos;
}

float Agent::getYPos()
{
    return y_pos;
}

float Agent::getZPos()
{
    return z_pos;
}

float Agent::getTime()
{
    return time;
}

float Agent::getSpeed()
{
    return speed;
}

float Agent::getPhiAngle()
{
    return phi_angle;
}

float Agent::getPsiAngle()
{
    return psi_angle;
}

void Agent::setThetaAngle( float angle)
{
    theta_angle = angle;
}

void Agent::setPhiAngle( float angle )
{
    phi_angle = angle;
}

void Agent::setPsiAngle( float angle )
{
    psi_angle = angle;
}

float Agent::getThetaAngle()
{
    return theta_angle;
}

void Agent::setSpeed( float s )
{
    speed = s;
}

void Agent::setID( int i )
{
    id = i;
}

void Agent::returnToStartingPosition()
{
    x_pos = starting_x;
    y_pos = starting_y;
    z_pos = starting_z;
    time = starting_time;
}

Coordinate* Agent::getCoordinates()
{
    return new Coordinate(x_pos, y_pos, z_pos, time);
}

float Agent::getStartingXPos(){ return starting_x; }
float Agent::getStartingYPos(){ return starting_y; }
float Agent::getStartingZPos(){ return starting_z; }
float Agent::getStartingTime(){ return starting_time; }


vector< pair<float,float> > Agent::getDisplacements()
{
    //cout << path.size() << endl;
    vector<pair<float, float> > displacements;
    int size = path.size();
    for ( int i = 0; i < size; i++ )
    {
        float diff_x = starting_x-path[i]->getX();
        float diff_y = starting_y-path[i]->getY();
        float diff_z = starting_z-path[i]->getZ();
        float time = path[i]->getTime();

        float dist = sqrt(diff_x*diff_x+diff_y*diff_y+diff_z*diff_z);
        displacements.push_back(make_pair(time, dist));

    }

    return displacements;
}

