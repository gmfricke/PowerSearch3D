#ifndef RECTANGULARPRISM_H
#define RECTANGULARPRISM_H

#include <boost/random/mersenne_twister.hpp>

#include "SearchSpace.h"
#include "rand_gen.h"
#include <vector>

class Coordinate;

class RectangularPrism : public SearchSpace
{
public:
    RectangularPrism( float xmin, float xmax, float ymin, float ymax, float zmin, float zmax, int nt);
    void updateSearcherPositions();
    float getVolume();

    float getHeight();
    float getDepth();
    float getWidth();

private:

};

#endif // RECTANGULARPRISM_H
