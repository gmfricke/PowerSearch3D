#ifndef COORDINATE_H
#define COORDINATE_H

class Coordinate
{
public:
    Coordinate();
    Coordinate(float x, float y, float z, float t);
    float getX();
    float getY();
    float getZ();
    float getTime();
    void setX(float x);
    void setY(float y);
    void setZ(float z);
    void setTime(float t);
    float getDistance( Coordinate* end );
    float getAngle(Coordinate *mid, Coordinate* end );
    float getSpeed( Coordinate* end );

private:
    float x;
    float y;
    float z;
    float t;
};

#endif // COORDINATE_H
