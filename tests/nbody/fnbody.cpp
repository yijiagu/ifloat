#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <immintrin.h>

#include <iostream>

#include <chrono>
using namespace std::chrono;

#define PI 3.141592653589793
#define SOLAR_MASS ( 4 * PI * PI )
#define DAYS_PER_YEAR 365.24

struct body {
   float x[3], fill, v[3], mass;
};

body solar_bodies[5];

static const int BODIES_SIZE = 5;

void offset_momentum(struct body *bodies, unsigned int nbodies)
{
   unsigned int i, k;

   
   for (i = 0; i < nbodies; ++i)
     for (k = 0; k < 3; ++k){

       bodies[0].v[k] = bodies[0].v[k] - bodies[i].v[k] * bodies[i].mass
	 / SOLAR_MASS;

     }
}

void bodies_advance(struct body *bodies, unsigned int nbodies, double dt)
{
   unsigned int N = (nbodies - 1) * nbodies / 2;
   static struct {
      float dx[3], fill;
   } r[1000];

   float mag[1000];
   unsigned int i, j, k, m;

   float dx[3][2], dsquared[2], distance[2], dmag[2];
     
   for(k = 0, i = 0; i < nbodies - 1; ++i)
      for(j = i + 1; j < nbodies; ++j, ++k)
         for ( m = 0; m < 3; ++m)
            r[k].dx[m] = bodies[i].x[m] - bodies[j].x[m];

   for (i = 0; i < N; i += 2) {
      for (m = 0; m < 3; ++m) {

         dx[m][0] = r[i].dx[m];

         dx[m][1] = r[i+1].dx[m];
      }



      dsquared[0] = dx[0][0] * dx[0][0] + dx[1][0] * dx[1][0] + dx[2][0] * dx[2][0];
      dsquared[1] = dx[0][1] * dx[0][1] + dx[1][1] * dx[1][1] + dx[2][1] * dx[2][1];

      distance[0] = sqrt(dsquared[0]);
      distance[1] = sqrt(dsquared[1]);


      dmag[0] = dt/(dsquared[0]*distance[0]);
      dmag[1] = dt/(dsquared[1]*distance[1]);
      
      mag[i] = dmag[0];
      mag[i+1] = dmag[1];
   }

   for (i = 0, k = 0; i < nbodies - 1; ++i)
      for ( j = i + 1; j < nbodies; ++j, ++k)
         for ( m = 0; m < 3; ++m) {

            bodies[i].v[m] = bodies[i].v[m] - r[k].dx[m] * bodies[j].mass
               * mag[k];

	    

            bodies[j].v[m] = bodies[j].v[m] + r[k].dx[m] * bodies[i].mass
               * mag[k];

         }

   for (i = 0; i < nbodies; ++i)
     for ( m = 0; m < 3; ++m){

        bodies[i].x[m] = bodies[i].x[m] + dt * bodies[i].v[m];

     }
}

float bodies_energy(struct body *bodies, unsigned int nbodies) {
   float dx[3], distance, e = 0.0;
   unsigned int i, j, k;

   for (i=0; i < nbodies; ++i) {
     float temp;
     

     temp = bodies[i].v[0] * bodies[i].v[0]
         + bodies[i].v[1] * bodies[i].v[1]
	+ bodies[i].v[2] * bodies[i].v[2];

       

     e = e + bodies[i].mass * temp * 0.5;

     
     for (j=i+1; j < nbodies; ++j) {
       for (k = 0; k < 3; ++k)
	 dx[k] = bodies[i].x[k] - bodies[j].x[k];

       float temp1;

       temp1 = dx[0] * dx[0] + dx[1] * dx[1] 
	 + dx[2] * dx[2];

	 
       distance = sqrt(temp1);
       float temp2;
       temp2 = bodies[i].mass * bodies[j].mass;
       e = e - temp2 / distance;
     }
   }
   return e;
}

int main()
{
  int i;
  int n = 1;
  
  solar_bodies[0].x[0] = 0.;
  solar_bodies[0].x[1] = 0.;
  solar_bodies[0].x[2] = 0.;
  solar_bodies[0].v[0] = 0.;
  solar_bodies[0].v[1] = 0.;
  solar_bodies[0].v[2] = 0.;
  solar_bodies[0].mass = SOLAR_MASS;

  solar_bodies[1].x[0] = 4.84143144246472090e+00;
  solar_bodies[1].x[1] = -1.16032004402742839e+00;
  solar_bodies[1].x[2] = -1.03622044471123109e-01; 
  solar_bodies[1].v[0] = 1.66007664274403694e-03 * DAYS_PER_YEAR;
  solar_bodies[1].v[1] = 7.69901118419740425e-03 * DAYS_PER_YEAR;
  solar_bodies[1].v[2] = -6.90460016972063023e-05 * DAYS_PER_YEAR;
  solar_bodies[1].mass = 9.54791938424326609e-04 * SOLAR_MASS;

  solar_bodies[2].x[0] = 8.34336671824457987e+00;
  solar_bodies[2].x[1] = 4.12479856412430479e+00;
  solar_bodies[2].x[2] = -4.03523417114321381e-01; 
  solar_bodies[2].v[0] = -2.76742510726862411e-03 * DAYS_PER_YEAR;
  solar_bodies[2].v[1] = 4.99852801234917238e-03 * DAYS_PER_YEAR;
  solar_bodies[2].v[2] = 2.30417297573763929e-05 * DAYS_PER_YEAR; 
  solar_bodies[2].mass = 2.85885980666130812e-04 * SOLAR_MASS;

  solar_bodies[3].x[0] =  1.28943695621391310e+01;
  solar_bodies[3].x[1] = -1.51111514016986312e+01;
  solar_bodies[3].x[2] = -2.23307578892655734e-01;
  solar_bodies[3].v[0] = 2.96460137564761618e-03 * DAYS_PER_YEAR;
  solar_bodies[3].v[1] = 2.37847173959480950e-03 * DAYS_PER_YEAR;
  solar_bodies[3].v[2] = -2.96589568540237556e-05 * DAYS_PER_YEAR;
  solar_bodies[3].mass = 4.36624404335156298e-05 * SOLAR_MASS;

  solar_bodies[4].x[0] = 1.53796971148509165e+01;
  solar_bodies[4].x[1] = -2.59193146099879641e+01;
  solar_bodies[4].x[2] =  1.79258772950371181e-01;
  solar_bodies[4].v[0] = 2.68067772490389322e-03 * DAYS_PER_YEAR;
  solar_bodies[4].v[1] = 1.62824170038242295e-03 * DAYS_PER_YEAR;
  solar_bodies[4].v[2] = -9.51592254519715870e-05 * DAYS_PER_YEAR;
  solar_bodies[4].mass = 5.15138902046611451e-05 * SOLAR_MASS;
  
  high_resolution_clock::time_point t_start = high_resolution_clock::now();
  for (i = 0; i < n; ++i)
    bodies_advance(solar_bodies, BODIES_SIZE, 0.015625);
  float energy;
  energy = bodies_energy(solar_bodies, BODIES_SIZE);
  
  high_resolution_clock::time_point t_end = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>( t_end - t_start).count();
  std::cout << "duration " << duration << std::endl; 

  std::cout << energy << std::endl;
  
  return 0;
}
