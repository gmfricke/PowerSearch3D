#include "SearchSpace.h"
#include "Target.h"
#include "Searcher.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdlib.h>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/lognormal_distribution.hpp>
#include <boost/math/distributions/pareto.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/random.hpp>
#include <vector>

using namespace std;
using namespace boost;
using namespace boost::math;
using boost::math::normal;
using namespace boost::numeric::ublas;

SearchSpace::SearchSpace( int nt )
{


    time_resolution = 1;
    n_searchers = 0;
    n_targets = nt;
    targets = new Target[n_targets];
    searchers = new Searcher[n_searchers];
    total_distance_travelled = 0;

    x_min_bound = -1;
    y_min_bound = -1;
    z_min_bound = -1;

    x_max_bound = 1;
    y_max_bound = 1;
    z_max_bound = 1;

    search_type = 0;
}

double SearchSpace::getTotalDistanceTravelled()
{
    float total_distance_travelled = 0;
    Searcher* searchers = getSearchers();
        for (int i = 0; i < getNumSearchers(); i++)
            total_distance_travelled += searchers->getDistanceTravelled();

  return total_distance_travelled;//*time_resolution;
}

double SearchSpace::getTotalTimeTravelled()
{
    float total_time_travelled = 0;
    Searcher* searchers = getSearchers();
        for (int i = 0; i < getNumSearchers(); i++)
            total_time_travelled += searchers->getTimeTravelled();

        //cout << "SearchSpace::getTotalTimeTravelled() time_resolution = " << time_resolution << endl;

        return total_time_travelled*time_resolution;
}

// Test whether any target is within detection radius of the line. Uses a commonly known point to line segment formula.
// start is the begining of the line segment, end is the end of the line segment we are interested in finding the distance to.
std::vector<Target*> SearchSpace::testpath( float start_x, float start_y, float start_z, float end_x, float end_y, float end_z )
{
    std::vector<Target*> found_targets;

  float x1 = start_x;
  float y1 = start_y;
  float z1 = start_z;

  float x2 = end_x;
  float y2 = end_y;
  float z2 = end_z;

  if ( norm(x2-x1,y2-y1,z2-z1) < 0.000001 ) return found_targets;



  // Iterate over targets and check distance from path vector
  for (int i = 0; i < n_targets; i++)
    {
      float d = numeric_limits<float>::max();

      float x0 = targets[i].getXPos();
      float y0 = targets[i].getYPos();
      float z0 = targets[i].getZPos();


      // Check if the closest approach is an endpoint of the line segment
      float vx = end_x - x1;
      float vy = end_y - y1;
      float vz = end_z - z1;

      float wx = x0 - x2;
      float wy = y0 - y2;
      float wz = z0 - z2;

      float c1 = wx*vx+wy*vy+wz*vz;
      float c2 = vx*vx+vy*vy+vz*vz;

      if ( c1  <= 0 )  // before segment start
      {
          d = norm(x0-x1, y0-y1, z0-z1);
      }
      else
      if ( c2 <= c1 ) // after segment end
      {
          d = norm(x0-x2, y0-y2, z0-z2);
      }
      else
      {
      // Else the closest approach is on the line segment

      float rx, ry, rz;
      cross_prod( x0-x1, y0-y1, z0-z1, x0-x2, y0-y2, z0-z2, rx, ry, rz );

      d = norm(rx, ry, rz)/norm(x2-x1,y2-y1,z2-z1);
        }
      //cout << "Detect radius: " << detection_radius << endl;
      if (d <= detection_radius)
	{
	  //cout << "Distance to point: " << d << ", detection radius: " << detection_radius << endl;
	  //cout << "Target Found" << endl;
	  targets[i].setFound();
      found_targets.push_back(targets+i);
	}
      
    }

  return found_targets;
}

void SearchSpace::cross_prod(float ux, float uy, float uz, float vx, float vy, float vz, float& rx, float& ry, float& rz)
{
  rx = uy*vz-uz*vy;
  ry = uz*vx-ux*vz;
  rz = ux*vy-uy*vx;
}

float SearchSpace::dot_prod(float x0, float y0, float z0, float x1, float y1, float z1)
{
  return x0*x1+y0*y1+z0*z1;
}

float SearchSpace::norm(float x, float y, float z)
{
  return sqrt(x*x+y*y+z*z);
}

Target* SearchSpace::getTargets()
{
    return targets;
}

Target* SearchSpace::getTargetsCopy(int& n)
{
    n = n_targets;
    Target* copy = new Target[n_targets];

    for (int i = 0; i < n_targets; i++) copy[i] = targets[i];
    return copy;
}

Searcher* SearchSpace::getSearchersCopy(int& n)
{
    n = n_searchers;
    Searcher* copy = new Searcher[n_searchers];
    for (int i = 0; i < n_searchers; i++) copy[i] = searchers[i];
    return copy;
}

Searcher* SearchSpace::getSearchers()
{
    return searchers;
}

int SearchSpace::getNumSearchers()
{
    return n_searchers;
}

int SearchSpace::getNumTargets()
{
    return n_targets;
}

void SearchSpace::setNumTargets(int v)
{
    n_targets = v;
}

SearchSpace::~SearchSpace()
{
  delete [] targets;
  delete [] searchers;
}

void SearchSpace::setSearchers( Searcher* s, int n )
{
    delete [] searchers;
    searchers = s;
    n_searchers = n;
}

void SearchSpace::setTargets( Target* s, int n )
{
    delete [] targets;
    targets = s;
    n_targets = n;
}

void SearchSpace::setBounds(float min_x, float min_y, float min_z, float max_x, float max_y, float max_z )
{
 x_min_bound = min_x;
 y_min_bound = min_y;
 z_min_bound = min_z;
 x_max_bound = max_x;
 y_max_bound = max_y;
 z_max_bound = max_z;

 // Fake 2d space
 if ( abs(min_x-max_x) < 0.00001 ) x_max_bound = min_x + 1;
 if ( abs(min_y-max_y) < 0.00001 ) y_max_bound = min_y + 1;
 if ( abs(min_z-max_z) < 0.00001 ) z_max_bound = min_z + 1;


}

float SearchSpace::getXMinBound()
{
    return x_min_bound;
}

float SearchSpace::getYMinBound()
{
    return y_min_bound;
}

float SearchSpace::getZMinBound()
{
    return z_min_bound;
}

float SearchSpace::getXMaxBound()
{
    return x_max_bound;
}

float SearchSpace::getYMaxBound()
{
    return y_max_bound;
}

float SearchSpace::getZMaxBound()
{
    return z_max_bound;
}

void SearchSpace::setTargetDetectRadius(float v)
{
    detection_radius = v;
}


float SearchSpace::getFromBetaPDF()
{
  float mean_scale = 6.9606; // So average movment size is comprable to brownian
  float mle_estParam = 1.2564e+6;
  double alpha = 1.843, beta = 1.5913e+7, randFromUnif;
  randFromUnif = rand()*1.0/RAND_MAX;
  beta_distribution<> dist(alpha, beta);
  double randFromDist = quantile(dist, randFromUnif);
  return randFromDist*mle_estParam*mean_scale;
}

float SearchSpace::getEuclideanDistance( Agent* a1, Agent* a2 )
{
    float x_delta = a1->getXPos()-a2->getXPos();
    float y_delta = a1->getYPos()-a2->getYPos();
    float z_delta = a1->getZPos()-a2->getZPos();
    return x_delta*x_delta+y_delta*y_delta+z_delta*z_delta;
}

float SearchSpace::getFromGenParetoPDF()
{
  float k = -0.12; //shape
  float sigma = 0.13; //scale
  float theta = 0.0006; //min value (location)
  float U; // Uniform random variable

  U = (float)rand()/(float)(RAND_MAX); //+1)?
  return theta + (sigma*(pow(U, k)-1))/(-k);
}

float SearchSpace::getFromExponentialPDF()
{
  float beta = 0.1033;
  float lambda = 1/beta;
  float U; // Uniform random variable
  U = (float)rand()/(float)(RAND_MAX); //+1)?
  return -log(U)/lambda;
}

float SearchSpace::getFromGammaPDF()
{
    return getFromGammaPDF(1.5914, 0.0730);
}

float SearchSpace::getFromGammaPDF(float shape, float scale)
{
  boost::gamma_distribution<> nd(shape, scale);

  boost::variate_generator<mt19937&, boost::gamma_distribution<> > var_gam(rng, nd);

  return var_gam();
  //return 1.0;
}

float SearchSpace::getFromNormalPDF(float mu, float sigma)
{
  boost::normal_distribution<> nd(mu, sigma);

  boost::variate_generator<mt19937&, boost::normal_distribution<> > var_nor(rng, nd);

  return var_nor();
}

float SearchSpace::getFromNormalPDF()
{
    return getFromNormalPDF(0.11, 0.0874);
}

float SearchSpace::getFromPowerLawPDF()
{
    // Mathamatica
    //pow((pow(x_min, alpha+1)-pow(x_max, alpha+1))*unif_var+pow(x_min, alpha+1),1.0/(alpha+1));

    float unif_var = rand()*1.0/RAND_MAX;
    float scale = 1;//2915;
    float exponent = 2.15;//2.908e-12;
    return pow(scale*unif_var, exponent);
}

float SearchSpace::getFromLogNormalPDF()
{
    float normal_var = getFromNormalPDF( -2.5, 0.95);
    return exp(normal_var);

    /*
  float sigma = 0.9570;
  float mu = 2.5;//-2.5;
  boost::lognormal_distribution<> nd(mu, sigma);

  boost::variate_generator<mt19937&, boost::lognormal_distribution<> > var_nor(rng, nd);

  return var_nor();
  */


}

float** SearchSpace::matmul4by4(float A[4][4], float B[4][4])
{
    float** C = new float*[4];
    for (int i = 0; i < 4; i++) C[i]=new float[4];

    for(int i = 0; i < 4; i++ )
        for(int j = 0; j < 4; j++ )
        {
            float rowcol_sum = 0;
            for (int k = 0; k < 4; k++ )
                rowcol_sum += A[i][k]*B[k][j];

            C[i][j] = rowcol_sum;
        }

    return C;
}

void SearchSpace::printMatrix(float** M)
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
            cout << M[i][j] << " ";
        cout << endl;
    }
    cout << endl;
}

void SearchSpace::clearTargets()
{
    n_targets = 0;
    delete [] targets;
    targets = new Target[0];
}

void SearchSpace::clearSearchers()
{
    n_searchers = 0;
    delete [] searchers;
    searchers = new Searcher[0];
}

void SearchSpace::setSearchType( int v )
{
    search_type = v;
}

int SearchSpace::getSearchType()
{
    return search_type;
}

string SearchSpace::getSearchTypeName()
{
    switch( search_type )
    {
    case 0:
        return "Observed";
        break;
    case 1:
        return "Brownian";
        break;
    case 2:
        return "LogNormal";
        break;
    case 3:
        return "GenPareto";
        break;
    case 4:
        return "PowerLaw";
        break;
    case 5:
        return "Exponential";
        break;
    case 6:
        return "CRW";
        break;
    case 7:
        return "CMLW";
        break;
    default:
        return "Unknown";
    }
}

void SearchSpace::setTimeResolution(float v)
{
    time_resolution = v;
}

float SearchSpace::getBoundingHeight()
{
    return abs( y_min_bound - y_max_bound );
}

float SearchSpace::getBoundingDepth()
{
    return abs( z_min_bound - z_max_bound );
}

float SearchSpace::getBoundingWidth()
{
    return abs( x_min_bound - x_max_bound );
}

float SearchSpace::angle(float a_x, float a_y, float a_z, float b_x, float b_y, float b_z)
{

    float a_norm = norm(a_x, a_y, a_z);
    float b_norm = norm(b_x, b_y, b_z);

    float unit_a_x = a_x/a_norm, unit_a_y = a_y/a_norm, unit_a_z = a_z/a_norm;
    float unit_b_x = b_x/b_norm, unit_b_y = b_y/b_norm, unit_b_z = b_z/b_norm;

    //cout << "Norms: " << norm(unit_a_x, unit_a_y, unit_a_z) << ", " << norm(unit_b_x, unit_b_y, unit_b_z) << endl;

    //return acos(dot_prod(unit_a_x, unit_a_y, unit_a_z, unit_b_x, unit_b_y, unit_b_z));

    float cross_x, cross_y, cross_z;
    cross_prod(unit_a_x, unit_a_y, unit_a_z, unit_b_x, unit_b_y, unit_b_z, cross_x, cross_y, cross_z);
    return atan2(norm(cross_x, cross_y, cross_z), dot_prod(unit_a_x, unit_a_y, unit_a_z,  unit_b_x, unit_b_y, unit_b_z));
}

// Rotate the input vector about the axis defined by the second input vector
// Uses: http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.html
Coordinate SearchSpace::rotateAboutVector(Coordinate vect, Coordinate axis, float angle)
{
    Coordinate result;

    float axis_L = norm(axis.getX(),axis.getY(),axis.getZ());
    Coordinate unit_axis;
    unit_axis.setX( axis.getX()/axis_L );
    unit_axis.setY( axis.getY()/axis_L );
    unit_axis.setZ( axis.getZ()/axis_L );

    float u = unit_axis.getX();
    float v = unit_axis.getY();
    float w = unit_axis.getZ();

    float x = vect.getX();
    float y = vect.getY();
    float z = vect.getZ();

    float cos_angle = cos(angle);
    float sin_angle = sin(angle);

    float dp = dot_prod(x,y,z,u,v,w);

    float x_row = u*dp*(1-cos_angle)+x*cos_angle+(-w*y+v*z)*sin_angle;
    float y_row = v*dp*(1-cos_angle)+y*cos_angle+(w*x-u*z)*sin_angle;
    float z_row = w*dp*(1-cos_angle)+z*cos_angle+(-v*x+u*y)*sin_angle;

    result.setX(x_row);
    result.setY(y_row);
    result.setZ(z_row);

    return result;
}

std::vector<double> SearchSpace::getFirstContactTimes()
{
    std::vector<double> first_contacts;
    for (int i = 0; i < n_searchers; i++)
    {
        int first_contact = searchers[i].getFirstTargetContactTime();
        if (first_contact > 1000000.0f-100.0f) continue; // Filter searchers who never found anything.
       // cout << "SearchSpace::getFirstContactTimes(): First contact: " << first_contact << endl;
        first_contacts.push_back( first_contact );
    }

 //   float min_time = INT_MAX;
 //  for (int i = 0; i < n_searchers; i++)
 //   {
 //       float time = searchers[i].getFirstTargetContactTime();
 //       if (min_time > time) min_time = time;
 //   }

    return first_contacts;
}

void SearchSpace::setCurrentTimeStep(int ts)
{
    current_time_step = ts;
}

int SearchSpace::getCurrentTimeStep()
{
    return current_time_step;
}

float SearchSpace::getCurrentTime()
{
    return getTimeResolution()*getCurrentTimeStep();
}

float SearchSpace::getTimeResolution()
{
    return time_resolution;
}

std::vector<double> SearchSpace::getSpeeds()
{
    std::vector<double> velocities;

    for ( int i = 0; i < n_searchers; i++ )
    {
        Searcher s = searchers[i];
        std::vector<Coordinate*> path = s.getPath();

        if (s.getPath().size() > 3)
        for (int j = 0; j < path.size()-2; j++)
        {
                Coordinate* start = path[j];
                Coordinate* middle = path[j+1];
                Coordinate* end = path[j+2];

                double step_length = norm(start->getX()-middle->getX(),start->getY()-middle->getY(),start->getZ()-middle->getZ());
                double velocity = step_length/time_resolution;
                velocities.push_back(velocity);
        }
    }

    return velocities;

}

std::vector<double> SearchSpace::getAngles()
{
    std::vector<double> angles;

    for ( int i = 0; i < n_searchers; i++ )
    {
        Searcher s = searchers[i];
        std::vector<Coordinate*> path = s.getPath();

        if (s.getPath().size() > 3)
        for (int j = 0; j < path.size()-2; j++)
        {

                Coordinate* start = path[j];
                Coordinate* middle = path[j+1];
                Coordinate* end = path[j+2];

                float a = angle(middle->getX() - start->getX(), middle->getY() - start->getY(), middle->getZ() - start->getZ(), end->getX() - middle->getX(), end->getY() - middle->getY(), end->getZ() - middle->getZ());
                angles.push_back(a);
        }
    }

return angles;
}

float SearchSpace::getSpeedFromObserved()
{
    // Choose uniformly from all observed speeds
    int n = observed_speeds.size();
    int choice = rand()%n;
    return observed_speeds[choice];
}



float SearchSpace::getAngleFromObserved()
{
    // Choose uniformly from all observed speeds
    int n = observed_angles.size();
    int choice = rand()%n;
    return observed_angles[choice];
}

void SearchSpace::setObservedAngles(std::vector<double> angles)
{
    observed_angles.clear();
    observed_angles = angles;
}

std::vector<double> SearchSpace::getObservedAngles()
{
    return observed_angles;
}

void SearchSpace::setObservedSpeeds(std::vector<double> speeds)
{
    observed_speeds.clear();
    observed_speeds = speeds;
}

std::vector<double> SearchSpace::getObservedSpeeds()
{
    return observed_speeds;
}






