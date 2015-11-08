#include "rand_gen.h"

#include <stdlib.h>

rand_gen::rand_gen(){

srand(time(NULL));
}
double rand_gen::uniform_dist(double a, double b){

double r=rand()/(RAND_MAX+1.0);
return a+(b-a)*r;
}

/*
To simulate exponential distribution exp(-lambda*x), the inverse
 method is used.
The cumulative distribution function for the exponential
distribution is:1-exp(-lambda*x). The inversion function
is  -log(1-U)/lambda. The simplified form is -log(U)/lambda,
where U is a uniform [0,1] random variable.
 http://cg.scs.carleton.ca/~luc/rnbookindex.html
*/
double rand_gen::exponential_dist(double lambda){
return(-1*log(uniform_dist(0.,1.))/lambda);
}

/*
Simulate power law with cut-off x^(-alpha)*exp(-lambda*x)
To simulate power-law with cutoff , one can generate an
exponentially distributed random number using the formula
above     (as k>0 and integer, so k start at 1)
and then accept or reject it with probability p or 1 - p
respectively (i.e accept U1 <p or reject U1>p, and
 U1 is a uniform [0,1] random variable),
where p = (x/x_min)^(-alpha)  and  x_min=1.

http://www.santafe.edu/~aaronc/powerlaws/
*/
double rand_gen::powerlaw_dist(double alpha, double lamda)
{

double x;
do{
x = exponential_dist(lamda);
} while (pow(x,-1*alpha) < uniform_dist(0.,1.));

return (x);

}
int rand_gen::binomial_dist(int n, double p){

int m=0;
for( int i=0;i<n;i++){

if(uniform_dist(0.,1.)<p) m++;
}
return m;
}

/*
Poisson distribution
f(n, lambda) = lambda^n exp(-lambda)/n!
Algorithm Poisson random number (Knuth):
init:          Let L =exp(-lambda), k=0 and p =1
do:          k = k+
Generate uniform random number u in [0,1]
and let p =p*u
while p > L.
return k - 1.
*/
int rand_gen::poisson_dist(double lambda){

double L=exp(-1*lambda);


int k;
double p=1.;
k = 0;

do{
k = k + 1;

p = p*uniform_dist(0.,1.);
} while (p > L);

return k - 1;
}
double rand_gen::normal_dist(double sigma, double mu)
{

double V1, V2, S;
double X;

do {
double U1 = uniform_dist(0.,1.);

double U2 = uniform_dist(0.,1.);

V1 = 2 * U1 - 1;

V2 = 2 * U2 - 1;
S = V1 * V1 + V2 * V2;
} while(S >= 1 || S == 0);

X = V1 * sqrt(-2 * log(S) / S);

return sigma*X+mu;
}
