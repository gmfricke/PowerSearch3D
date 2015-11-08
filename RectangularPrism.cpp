#include "RectangularPrism.h"
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/lognormal_distribution.hpp>
#include <boost/math/distributions/pareto.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/random.hpp>
#include "Searcher.h"
#include "Coordinate.h"
#include "Agent.h"

#define PI 3.141592653589793

using namespace boost;
using namespace boost::math;
using boost::math::normal;
using namespace boost::numeric::ublas;

using namespace std;

RectangularPrism::RectangularPrism( float xmin, float xmax, float ymin, float ymax, float zmin, float zmax, int nt ) : SearchSpace(nt)
{
    x_min_bound = xmin;
    x_max_bound = xmax;
    y_min_bound = ymin;
    y_max_bound = ymax;
    z_min_bound = zmin;
    z_max_bound = zmax;
}

void RectangularPrism::updateSearcherPositions()
{
    float min_step = 0.0001; // to avoid flaoting point zero errors
    float scale = 0.0;
    float angle = 0.0;
    bool corralate_turning_angles = false;
    bool bootstrap_angles = false;
    float change_in_inclination;
    float change_in_azmuth;
    float old_inclination;
    float old_azmuth;
    float old_radius;
    Coordinate new_vect, previous_vector;

        //cout << "Moving " << n_searchers << endl;


        for (int i = 0; i < n_searchers; i++ )
        {
            switch( search_type )
            {
            case 0:
                cout << "RectangularPrism: This function should not be called when the search type is set to \"Observed\"" << endl;
                return;
                break;
            case 1: scale = 0.11*time_resolution;//getFromMaxwellPDF()*time_resolution;
                        //getFromNormalPDF()*time_resolution;
                //cout << "RectangularPrism:updateSearcherPositions: Brownian" << endl;
                break;
            case 2: scale = getFromLogNormalPDF();
                //cout << scale << endl;
                scale *= time_resolution;
                //cout << "RectangularPrism:updateSearcherPositions: LogNormal" << endl;
                break;
            case 3: scale = getFromGenParetoPDF()*time_resolution;
               // scale *=  6.9606; // So mean velocity is comparable to Brownian
                //cout << "RectangularPrism:updateSearcherPositions: Gen. Pareto" << endl;
                break;
            case 4: scale = getFromPowerLawPDF()*time_resolution;//1.0, 70, -2.15)*time_resolution;
                //cout << "RectangularPrism:updateSearcherPositions: Power Law, scale = " << scale << endl;
                //scale *= 0.0079;//scale *= 0.4104; //1.64179;//0.021; // Make mean equal to observed mean speed of 0.11
                //cout << scale << endl;
                break;
            case 5: scale = getFromExponentialPDF()*time_resolution;
                //cout << "RectangularPrism:updateSearcherPositions: Exponential, scale = " << scale << endl;
                break;
            case 6:
                scale = getFromGammaPDF(1.5914, 0.0730);
                //cout << "RectangularPrism:updateSearcherPositions: Gamma, scale = " << scale << endl;
                scale *= time_resolution;
                break;
            case 7: scale = 0.11*time_resolution; // 0.11 is the WT mean velocity
               // cout << "RectangularPrism:updateSearcherPositions: Corralated Random Walk, scale = " << scale << endl;
                corralate_turning_angles = true;
                break;
            case 8: scale = getFromLogNormalPDF()*time_resolution;
               // cout << "RectangularPrism:updateSearcherPositions: Corralation Modulated LogNormal Walk, scale = " << scale << endl;
                corralate_turning_angles = true;
                break;
            case 9: scale = getFromGammaPDF(1.5914, 0.0730)*time_resolution;
               // cout << "RectangularPrism:updateSearcherPositions: Corralation Modulated LogNormal Walk, scale = " << scale << endl;
                corralate_turning_angles = true;
                break;
            case 10:
                scale = getSpeedFromObserved()*time_resolution; // Bootstrap
                //angle = getAngleFromObserved();
                //cout << "RectangularPrism:updateSearcherPositions: Bootstrap, scale = " <<
                break;
            case 11:
                scale = getSpeedFromObserved()*time_resolution; // Bootstrap
                bootstrap_angles = true;
                break;

            default:
                cout << "RectangularPrism:updateSearcherPositions: unknown search type " << search_type << endl;
            }

            //cout << scale << endl;

            if (scale < min_step) scale = min_step;

            float x_pos = searchers[i].getXPos();
            float y_pos = searchers[i].getYPos();
            float z_pos = searchers[i].getZPos();

            float old_x = x_pos;
            float old_y = y_pos;
            float old_z = z_pos;


                // Results in a normal distribution of turning angles not uniform, but OK for first step since there is no previous step to correlate with
                //new_step_x = ((rand()%10000)*2.0/10000)-1.0;
                //new_step_y = ((rand()%10000)*2.0/10000)-1.0;
                //new_step_z = ((rand()%10000)*2.0/10000)-1.0;

            // Useful theorem
            //For the direction, a useful fact is that if V=<X1,X2,X3> where X1,X2,X3 are independent normal random variables with mean 0 and variance 1, then
            // V/norm(V) is uniformly distributed on the unit sphere (Muller 1959, Marsaglia 1972) Muller, M. E. "A Note on a Method for Generating Points Uniformly on N-Dimensional Spheres." Comm. Assoc. Comput. Mach. 2, 19-20, Apr. 1959.

            float new_step_x, new_step_y, new_step_z;
            std::vector<Coordinate*> path = searchers[i].getPath();
            int size = path.size();
            if (size < 2)
            {
                new_step_x = getFromNormalPDF(0,1); // since the std of var = 1 is sqrt(1) = 1
                new_step_y = getFromNormalPDF(0,1);
                new_step_z = getFromNormalPDF(0,1);
            }
            else
            {
                // move previous vector to origin

                float previous_vector_x = path[size-1]->getX() - path[size-2]->getX();
                float previous_vector_y = path[size-1]->getY() - path[size-2]->getY();
                float previous_vector_z = path[size-1]->getZ() - path[size-2]->getZ();
                //cout << "Previous vector: <" << previous_vector_x << "," << previous_vector_y << "," << previous_vector_z << ">" << endl;

                // Convert to spherical coordinates

                float radius = norm(previous_vector_x, previous_vector_y, previous_vector_z);
                float inclination = acos(previous_vector_z/radius);
                float azmuth = atan2(previous_vector_y, previous_vector_x);

                 old_radius = radius;
                 old_inclination = inclination;
                 old_azmuth = azmuth;

                //cout << "Previous vector: <" << previous_vector_x << "," << previous_vector_y << "," << previous_vector_z << ">" << endl;

               // cout << "Azmuth:" << azmuth << endl;
               // cout << "Inclination:" << inclination << endl;
               // cout << "Radius:" << radius << endl;



               if (bootstrap_angles)
               {
                   change_in_inclination = getAngleFromObserved();
               }
              else if (corralate_turning_angles)
               {
                   // Estimated from observation Gamma distribution with
                   //a = 2.1586 (shape)
                   //b = 34.3086 (scale)
                   change_in_inclination = getFromGammaPDF(2.1586,0.598799);//*PI; Don't need this "*PI" because the value is already in radians
                    //change_in_inclination = getFromNormalPDF(0,0.1)*PI;
                   // if (change_in_inclination > PI)
                   // {
                   //     change_in_inclination = 0;
                   // }
                }
               else
               {
                   change_in_inclination = (rand()*1.0/RAND_MAX)*PI;
               }

                change_in_azmuth = (rand()*1.0/RAND_MAX)*2*PI;

                inclination += change_in_inclination;
                //azmuth += change_in_azmuth;

                //new_step_x = previous_vector_x;
                //new_step_y = previous_vector_y;
                //new_step_z = previous_vector_z;

                // Convert to Cartesian coordinates
                new_step_x = radius*sin(inclination)*cos(azmuth);
                new_step_y = radius*sin(inclination)*sin(azmuth);
                new_step_z = radius*cos(inclination);

                // rotate about the z-axis (doesn't matter which axis we choose)

                new_vect.setX(new_step_x);
                new_vect.setY(new_step_y);
                new_vect.setZ(new_step_z);

                // cout << "Vector: <" << vect.getX() << "," << vect.getY() << "," << vect.getZ() << ">" << endl;

                previous_vector.setX(previous_vector_x);
                previous_vector.setY(previous_vector_y);
                previous_vector.setZ(previous_vector_z);

          //      Coordinate z_axis;
          //      z_axis.setX(0);
          //      z_axis.setY(0);
          //      z_axis.setZ(1);


                //vect = rotateAboutVector(vect, z_axis, change_in_inclination); // inclination for lack of a better name

            //    vect.setX( vect.getX()*cos(change_in_azmuth) - vect.getY()*sin(change_in_azmuth) );
            //    vect.setX( vect.getX()*sin(change_in_azmuth) + vect.getY()*cos(change_in_azmuth) );
                // z is unchanged

                //cout << "Rotated vector about z-axis by : " << change_in_inclination/PI << "\u03C0 <" << vect.getX() << "," << vect.getY() << "," << vect.getZ() << ">" << endl;


                // Now rotate the new vector using the previous vector as the axis
                new_vect = rotateAboutVector(new_vect, previous_vector, change_in_azmuth); // azmuth for lack of a better name

                //cout << "Rotated vector about previous by : " << change_in_azmuth/PI << "\u03C0 <" << vect.getX() << "," << vect.getY() << "," << vect.getZ() << ">" << endl;

                new_step_x = new_vect.getX();
                new_step_y = new_vect.getY();
                new_step_z = new_vect.getZ();
            }


            //cout << "<" << old_x << "," << old_y << ","<< old_z << "> <" << new_step_x << ","<< new_step_y << ","<< new_step_z << ">" << endl;
            float scale_step = norm(new_step_x, new_step_y, new_step_z);
            new_step_x /= scale_step;
            new_step_y /= scale_step;
            new_step_z /= scale_step;

            new_step_x = scale*new_step_x; // dividing by the scale_step makes the distribution uniform because normal distributions are radially symmetric and also makes it a unit vector
            new_step_y = scale*new_step_y; // Multiplyng by the scale gives the length drawn from one of the probability distributions above
            new_step_z = scale*new_step_z;

            x_pos = x_pos + new_step_x;
            y_pos = y_pos + new_step_y;
            z_pos = z_pos + new_step_z;

            //cout << "New vector: <" << x_pos << "," << y_pos << "," << z_pos << ">" << endl;
           // cout << angle(old_x, old_y, old_z, new_step_x, new_step_y, new_step_z) << endl;

            if ( !(x_pos < x_min_bound || x_pos > x_max_bound || y_pos < y_min_bound || y_pos > y_max_bound || z_pos < z_min_bound || z_pos > z_max_bound) )
            {

                /*
                float new_radius = norm(x_pos, y_pos, z_pos);
                float new_inclination = acos(z_pos/new_radius);
                float new_azmuth = atan2(y_pos, x_pos);

              cout << "------------" << endl;
              cout << "<" << old_x << "," << old_y << "," << old_z << ">, <" << x_pos << "," << y_pos << "," << z_pos << ">: ";
              cout << (180.0/PI)*angle(previous_vector.getX(), previous_vector.getY(), previous_vector.getZ(), new_vect.getX(), new_vect.getY(), new_vect.getZ()) << endl;
              cout << "Old inclination: " << (180.0/PI)*old_inclination << endl;
              cout << "Old azmuth: " << (180.0/PI)*old_azmuth << endl;
              cout << "Change in inclination: " << (180.0/PI)*change_in_inclination << endl;
              cout << "Change in azmuth: " << (180.0/PI)*change_in_azmuth << endl;
              cout << "New inclination: " << (180.0/PI)*new_inclination << endl;
              cout << "New azmuth: " << (180.0/PI)*new_azmuth << endl;
              cout << "------------" << endl;
              */

              searchers[i].setXPos(x_pos);
              searchers[i].setYPos(y_pos);
              searchers[i].setZPos(z_pos);
              searchers[i].setTime(getCurrentTime());

              if (boost::math::isnan(x_pos))
              {
                  cout << "RectangularPrism::updateSearcherPositions(): x_pos is NaN" << endl;
              }

              searchers[i].addToPath(new Coordinate(searchers[i].getXPos(),searchers[i].getYPos(),searchers[i].getZPos(),searchers[i].getTime())); // Do it here so we avoid concurrent access problems

              //cout << (180.0/PI)*angle(previous_vector.getX(), previous_vector.getY(), previous_vector.getZ(), new_vect.getX(), new_vect.getY(), new_vect.getZ()) << endl;

            }
    }

}

float RectangularPrism::getVolume()
{
    float width = x_max_bound - x_min_bound;
    float height = y_max_bound - y_min_bound;
    float depth = z_max_bound - z_min_bound;

    return width*depth*height;
}

float RectangularPrism::getHeight()
{
    return getBoundingHeight();
}

float RectangularPrism::getDepth()
{
    return getBoundingDepth();
}

float RectangularPrism::getWidth()
{
    return getBoundingWidth();
}
