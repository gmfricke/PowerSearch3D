#ifndef AGENT_H
#define AGENT_H

#include <vector>
#include <utility> // for pair class

using namespace std;

class Coordinate;

class Agent
{
public:
    Agent();
    void setXPos( float x );
    void setYPos( float y );
    void setZPos( float z );
    void setTime( float t );
    void setStartingXPos( float x );
    void setStartingYPos( float y );
    void setStartingZPos( float z );
    void setStartingTime( float t );

    float getStartingXPos();
    float getStartingYPos();
    float getStartingZPos();
    float getStartingTime();

    float getXPos();
    float getYPos();
    float getZPos();
    float getTime();
    float getPsiAngle();
    float getPhiAngle();
    float getThetaAngle();
    void setPsiAngle(float angle);
    void setPhiAngle(float angle);
    void setThetaAngle(float angle);
    float getRadius();
    float getSpeed();
    void setSpeed( float s );
    void setRadius(float r);
    int getID();
    void setID( int i );
    void returnToStartingPosition();
    Coordinate* getCoordinates();

    vector<pair<float, float> > getDisplacements(); // All displacements

    int red;
    int blue;
    int green;

protected:
    float starting_x;
    float starting_y;
    float starting_z;
    float starting_time;
    float x_pos;
    float y_pos;
    float z_pos;
    float time;
    float speed; // units?
    float psi_angle; // in radians
    float phi_angle;
    float theta_angle;
    float radius;
    int id;
    vector<Coordinate*> path;
};

#endif // AGENT_H
