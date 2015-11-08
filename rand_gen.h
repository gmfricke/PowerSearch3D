#ifndef RAND_GEN_H
#define RAND_GEN_H

#include<ctime>

#include<cmath>

using namespace std;


class rand_gen{
public:

rand_gen();
double uniform_dist(double a, double b);

double exponential_dist(double lamda);
double powerlaw_dist(double alpha, double lamda);

int binomial_dist(int n, double p);
int poisson_dist(double lambda);

double normal_dist(double sigma, double mu);
};

#endif // RAND_GEN_H
