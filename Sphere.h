#ifndef SEARCHSPHERE_H
#define SEARCHSPHERE_H

#include <boost/random/mersenne_twister.hpp>

#include "SearchSpace.h"
#include "rand_gen.h"

using namespace boost;

class Agent;

class Sphere : public SearchSpace
{
public:
    Sphere( float r, int nt, float dr, int sp );
    void updateTargetPositions();
    void updateSearcherPositionsBrownian();
    void updateSearcherPositionsLogNormal();
    void updateSearcherPositionsLevy();
    void updateSearcherPositionsBeta();
    void updateSearcherPositions();
    float getRadius();

    ~Sphere();

private:
    mt19937 rng;
    float radius;
    rand_gen generator;

};

#endif // SEARCHSPHERE_H
