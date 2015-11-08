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
#include "Sphere.h"
#include "Searcher.h"
#include "Coordinate.h"
#include "Agent.h"

using namespace boost;
using namespace boost::math;
using boost::math::normal;
using namespace boost::numeric::ublas;
using namespace std;

#define PI 3.14159

Sphere::Sphere( float r, int nt, float dr, int sp) : SearchSpace(nt)
{
  detection_radius = dr;
  radius = r;
  search_pattern = sp;
}


void Sphere::updateTargetPositions()
{

}

void Sphere::updateSearcherPositionsLogNormal()
{

    for (int i = 0; i < n_searchers; i++ )
    {
        float x_pos = searchers[i].getXPos();
        float y_pos = searchers[i].getYPos();
        float z_pos = searchers[i].getZPos();

        float new_step_x = ((rand()%10000)*2.0/10000)-1.0;
        float new_step_y = ((rand()%10000)*2.0/10000)-1.0;
        float new_step_z = ((rand()%10000)*2.0/10000)-1.0;

        float scale = getFromLogNormalPDF();

        float scale_step = sqrt(new_step_x*new_step_x+new_step_y*new_step_y+new_step_z*new_step_z)*100;

        new_step_x = scale*new_step_x/scale_step;
        new_step_y = scale*new_step_y/scale_step;
        new_step_z = scale*new_step_z/scale_step;

        x_pos = x_pos + new_step_x;
        y_pos = y_pos + new_step_y;
        z_pos = z_pos + new_step_z;

        if (sqrt(x_pos*x_pos + y_pos*y_pos +z_pos*z_pos) < radius)
        {
            searchers[i].setXPos(x_pos);
            searchers[i].setYPos(y_pos);
            searchers[i].setZPos(z_pos);
        }

    }
}


void Sphere::updateSearcherPositionsBeta()
{

    for (int i = 0; i < n_searchers; i++ )
    {
        float x_pos = searchers[i].getXPos();
        float y_pos = searchers[i].getYPos();
        float z_pos = searchers[i].getZPos();

        float new_step_x = ((rand()%10000)*2.0/10000)-1.0;
        float new_step_y = ((rand()%10000)*2.0/10000)-1.0;
        float new_step_z = ((rand()%10000)*2.0/10000)-1.0;

        float scale = getFromBetaPDF();

        float scale_step = sqrt(new_step_x*new_step_x+new_step_y*new_step_y+new_step_z*new_step_z)*100;

        new_step_x = scale*new_step_x/scale_step;
        new_step_y = scale*new_step_y/scale_step;
        new_step_z = scale*new_step_z/scale_step;

        x_pos = x_pos + new_step_x;
        y_pos = y_pos + new_step_y;
        z_pos = z_pos + new_step_z;

        if (sqrt(x_pos*x_pos + y_pos*y_pos +z_pos*z_pos) < radius)
        {
            searchers[i].setXPos(x_pos);
            searchers[i].setYPos(y_pos);
            searchers[i].setZPos(z_pos);
        }

    }
}

void Sphere::updateSearcherPositions()
{
    //cout << "Number of searchers: " << n_searchers << endl;

    for (int i = 0; i < n_searchers; i++ )
    {
        float x_pos = searchers[i].getXPos();
        float y_pos = searchers[i].getYPos();
        float z_pos = searchers[i].getZPos();

    float old_x = x_pos;
    float old_y = y_pos;
    float old_z = z_pos;

        float new_step_x = ((rand()%10000)*2.0/10000)-1.0;
        float new_step_y = ((rand()%10000)*2.0/10000)-1.0;
        float new_step_z = ((rand()%10000)*2.0/10000)-1.0;

    float scale = 0.0;

    switch(search_pattern)
      {
      case 1:
        scale = getFromBetaPDF();
        break;
      case 2:
        scale = getFromGammaPDF();
        break;
      case 3:
        scale = getFromLogNormalPDF();
        break;
      case 4:
        //scale = 0.1033;
        scale = getFromNormalPDF();
        break;
      case 5:
        scale = getFromGenParetoPDF();
        break;
      case 6:
        scale = getFromExponentialPDF();
        break;
      default:
        cout << "Sphere::updateSearcherPositions(): Invalid search pattern type: " << search_pattern << endl;
        exit(1);
        break;
      }
    //cout << scale << endl;

        float scale_step = sqrt(new_step_x*new_step_x+new_step_y*new_step_y+new_step_z*new_step_z)/100.0;

        new_step_x = scale*new_step_x/scale_step;
        new_step_y = scale*new_step_y/scale_step;
        new_step_z = scale*new_step_z/scale_step;

        x_pos = x_pos + new_step_x;
        y_pos = y_pos + new_step_y;
        z_pos = z_pos + new_step_z;

        //if (sqrt(x_pos*x_pos + y_pos*y_pos +z_pos*z_pos) < radius)
        {

      // Sets the found flag on targets within distance d of the searcher's path
      testpath(old_x, old_y, old_z, x_pos, y_pos, z_pos);

      total_distance_travelled += norm(old_x-x_pos, old_y-y_pos, old_z-z_pos);

          searchers[i].setXPos(x_pos);
          searchers[i].setYPos(y_pos);
          searchers[i].setZPos(z_pos);

          cout << "Step length: " << norm(old_x-x_pos, old_y-y_pos, old_z-z_pos) << endl;

          searchers[i].addToPath(new Coordinate(searchers[i].getXPos(),searchers[i].getYPos(),searchers[i].getZPos(),searchers[i].getTime())); // Do it here so we avoid concurrent access problems

        }

}
}






void Sphere::updateSearcherPositionsBrownian()
{
    for (int i = 0; i < n_searchers; i++ )
    {
        float x_pos = searchers[i].getXPos();
        float y_pos = searchers[i].getYPos();
        float z_pos = searchers[i].getZPos();
        //float phi = searchers[i].getPhiAngle();
        //float psi = searchers[i].getPsiAngle();
        //float theta = searchers[i].getThetaAngle();

        //cout << "(" << x_pos << ", " << y_pos << ", " << z_pos << ")";


        float speed = searchers[i].getSpeed();
        //float direction = searchers[i].getDirection();

        //float phi_delta = (float)(2*PI*rand()/(1.0f*RAND_MAX));
        //float psi_delta = (float)(2*PI*rand()/(1.0f*RAND_MAX));
        //float theta_delta = (float)(2*PI*rand()/(1.0f*RAND_MAX));

        //cout << "phi delta" << phi_delta << endl;
        //cout << "theta delta" << theta_delta << endl;
        //cout << "psi delta" << psi_delta << endl;


        //direction=direction+direction_delta;

        // OpenGL uses right handed rotation

        //float x_rotation[4][4] = {{1, 0, 0, 0},
        //                          {0, cos(phi_delta), sin(phi_delta), 0},
        //                          {0, -sin(phi_delta), cos(phi_delta), 0},
    //                           {0, 0, 0, 1}};


        //float y_rotation[4][4] = {{cos(theta_delta), 0, -sin(theta_delta), 0},
        //                          {0, 1, 0, 0},
        //                          {sin(theta_delta), 0, cos(theta_delta), 0},
                    //                           {0, 0, 0, 1}};

      //float z_rotation[4][4] = {{cos(psi_delta), sin(psi_delta), 0, 0},
      //                         {-sin(psi_delta), cos(psi_delta), 0, 0},
          //                        {0, 0, 1, 0},
      //                         {0, 0, 0, 1}};

        // Multiply the matricies
        //float **C, **total_rotation;
        //float total_rotation[4][4] = {{1,2,3,4},{5,6,7,8},{7,3,2,7},{9,8,7,6}};
        //C = matmul4by4(x_rotation, y_rotation);

        // Make a copy of the matrix
        //float tempM[4][4];
        //for (int k = 0; k < 4; k++ )
        //    for (int j = 0; j < 4; j++ )
        //       tempM[k][j] = C[k][j];

        //total_rotation = matmul4by4(tempM, z_rotation);

        //float position_vector[4] = {x_pos, y_pos, z_pos, 1};

        // Apply the rotation transformation to the position vector for this searcher
        //float new_vector[4];
        //new_vector[0] = (total_rotation[0][0]*position_vector[0]+total_rotation[0][1]*position_vector[1]+total_rotation[0][2]*position_vector[2]+total_rotation[0][3]*position_vector[3]);
        //new_vector[1] = (total_rotation[1][0]*position_vector[0]+total_rotation[1][1]*position_vector[1]+total_rotation[1][2]*position_vector[2]+total_rotation[1][3]*position_vector[3]);
        //new_vector[2] = (total_rotation[2][0]*position_vector[0]+total_rotation[2][1]*position_vector[1]+total_rotation[2][2]*position_vector[2]+total_rotation[2][3]*position_vector[3]);
        //new_vector[3] = (total_rotation[3][0]*position_vector[0]+total_rotation[3][1]*position_vector[1]+total_rotation[3][2]*position_vector[2]+total_rotation[3][3]*position_vector[3]);

        //x_pos = x_pos + speed * cos( psi );
        //y_pos = y_pos - speed * sin( psi );

        //float new_x_pos = new_vector[0];
        //float new_y_pos = new_vector[1];
        //float new_z_pos = new_vector[2];

        float new_step_x = ((rand()%10000)*2.0/10000)-1.0;
        float new_step_y = ((rand()%10000)*2.0/10000)-1.0;
        float new_step_z = ((rand()%10000)*2.0/10000)-1.0;

        //cout << new_step_x << ", " << new_step_y << ", " << new_step_z << endl;

        float scale_step = sqrt(new_step_x*new_step_x+new_step_y*new_step_y+new_step_z*new_step_z)*100;

        // Make the vector a unit vector
        new_step_x = new_step_x/scale_step;
        new_step_y = new_step_y/scale_step;
        new_step_z = new_step_z/scale_step;

        //float step_length = sqrt(new_step_x*new_step_x + new_step_y*new_step_y + new_step_z*new_step_z);

        //float new_length = sqrt(new_x_pos*new_x_pos + new_y_pos*new_y_pos +new_z_pos*new_z_pos);

        //float new_step_x = new_x_pos*speed;
        //float new_step_y = new_y_pos*speed;
        //float new_step_z = new_z_pos*speed;

        x_pos = x_pos + new_step_x;
        y_pos = y_pos + new_step_y;
        z_pos = z_pos + new_step_z;

//        cout << "Step scale: " << scale_step << endl;
//        cout << "Step length: " << step_length << endl;


       //cout << "Rotated: (" << x_pos << ", " << y_pos << ", " << z_pos << ")" << endl;

        //cout << " Coords: (" << x_pos << ", " << y_pos << ", " << z_pos << ")" << endl;

        //cout << length << endl;

        if (sqrt(x_pos*x_pos + y_pos*y_pos +z_pos*z_pos) < radius)
        {
        searchers[i].setSpeed(speed);
        searchers[i].setXPos(x_pos);
        searchers[i].setYPos(y_pos);
        searchers[i].setZPos(z_pos);
        }
        else
        {
            //cout << "Searcher " << i << ": " << sqrt(x_pos*x_pos + y_pos*y_pos) << "exceeded boundary " << radius << endl;
        }
    }
}


float Sphere::getRadius()
{
    return radius;
}

void placeTargets()
{

}

Sphere::~Sphere()
{

}
