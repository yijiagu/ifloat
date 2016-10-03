/*
*   GAUSSIAN TRIPLE INTEGRAL ALGORITHM 4.6
* 
*   To approximate I = triple integral ( ( f(x,y,z) dz dy dx ) ) with limits
*   of integration from a to b for x, from c(x) to d(x) for y, and from
*   alpha(x,y) to beta(x,y) for z.
*
*   INPUT:   endpoints a, b; postive integers m, n, p. (Assume that the
*            roots r(i,j) and coefficients c(i,j) are available for i
*            equals m, n, and p and for 1 <= j <= i.
*
*   OUTPUT:  approximation J TO I. 
*/

#include<stdio.h>
#include<math.h>
#include <iostream>

#include <InstableFloat.h>
#include <InstableEnv.h>
#include <InstableMath.h>

#include <chrono>
using namespace std::chrono;

#define true 1
#define false 0

InstableFloat F(InstableFloat, InstableFloat);
InstableFloat C(); 
InstableFloat D(InstableFloat); 
InstableFloat ALPHA(InstableFloat, InstableFloat);
InstableFloat BETA();

void INPUT(double *, double *);

int main(){
   double r[5][5],co[5][5];
   double A_float[10], B_float[10];
   
   InstableFloat A,B,H1,H2,AJ,JX,D1,C1,K1,K2,JY,Z1,Z2,L1,L2,X,Y,Z,Q; 
   int N,M,I,J,P,K,OK; 
   
   OK = 1;
   N = 5;
   M = 2;
   P = 5;
   
   INPUT(A_float, B_float);

   for(int i = 0; i < 1; i ++){

     high_resolution_clock::time_point t_start = high_resolution_clock::now();
     
     A = A_float[i];
     B = B_float[i];
     
     if (OK) {
       r[1][0] = 0.5773502692; r[1][1] = -r[1][0]; co[1][0] = 1.0;
       co[1][1] = 1.0; r[2][0] = 0.7745966692; r[2][1] = 0.0;
       r[2][2] = -r[2][0]; co[2][0] = 0.5555555556; co[2][1] = 0.8888888889;
       co[2][2] = co[2][0]; r[3][0] = 0.8611363116; r[3][1] = 0.3399810436;
       r[3][2] = -r[3][1]; r[3][3] = -r[3][0]; co[3][0] = 0.3478548451;
       co[3][1] = 0.6521451549; co[3][2] = co[3][1]; co[3][3] = co[3][0];
       r[4][0] = 0.9061798459; r[4][1] = 0.5384693101; r[4][2] = 0.0;
       r[4][3] = -r[4][1]; r[4][4] = -r[4][0]; co[4][0] = 0.2369268850;
       co[4][1] = 0.4786286705; co[4][2] = 0.5688888889; co[4][3] = co[4][1];
       co[4][4] = co[4][0];
       /* STEP 1 */
       InstableFloat t1, t2;
       t1 = B - A;
       t2 = B + A;
       
       H1 = t1 / 2.0;       
       H2 = t2 / 2.0;
       
       AJ = 0.0;                         /* Use AJ instead of J. */ 
       /* STEP 2 */ 
       for (I=1; I<=M; I++) {    
	 /* STEP 3 */
	 InstableEnv::open_label(1);
	 X = H1 * r[M-1][I-1] + H2;
	 InstableEnv::close_label();
	 
	 JX = 0.0;
	 C1 = C();
	 D1 = D(X);
	 
	 InstableFloat t5, t6;
	 t5 = D1 - C1;
	 t6 = D1 + C1;
	 K1 = t5 / 2.0;
	 K2 = t6 / 2.0;
	 
	 /* STEP 4 */ 
	 for (J=1; J<=N; J++) { 
	   /* STEP 5 */
	   InstableEnv::open_label(2);
	   Y = K1 * r[N-1][J-1] + K2;
	   InstableEnv::close_label();
	   
	   JY = 0.0;
	   /* use Z1 for Beta and Z2 for Alpha */
	   Z1 = BETA();
	   Z2 = ALPHA(X, Y);
	   InstableFloat t3, t4;
	   t3 = Z1 - Z2;
	   t4 = Z1 + Z2;
	   L1 = t3 / 2.0;
	   L2 = t4 / 2.0;
	   /* STEP 6 */ 
	   for (K=1; K<=P; K++) {
	     InstableEnv::open_label(3);
	     Z = L1 * r[P-1][K-1] + L2;
	     InstableEnv::close_label();
	     
	     Q = F(X, Y);
	     
	     InstableEnv::open_label(4);
	     JY = JY + co[P-1][K-1] * Q;
	     InstableEnv::close_label();
	   }
	   /* STEP 7 */
	   InstableEnv::open_label(5);
	   JX = JX + co[N-1][J-1] * L1 * JY;
	   InstableEnv::close_label();
	 }
	 /* STEP 8 */
	 InstableEnv::open_label(6);
	 AJ = AJ + co[M-1][I-1] * K1 * JX;
	 InstableEnv::close_label();
       }
       /* STEP 9 */ 
       AJ = AJ * H1;
       
       high_resolution_clock::time_point t_end = high_resolution_clock::now();
       auto duration = duration_cast<microseconds>( t_end - t_start).count();
       std::cout << "duration " << duration << std::endl; 
       /* STEP 10 */ 
       std::cout << AJ << std::endl;
     }
   } 
   return 0;
}

/* Change functions F,C,D,ALPHA,BETA for a new problem */
InstableFloat F(InstableFloat X, InstableFloat Y)
{  /* F is the integrand */ 
  InstableFloat f, t;
  InstableEnv::open_label(7);
  t = X * X + Y * Y;
  InstableEnv::close_label();
  f = t;
  
  //f = sqrt(t); 
  return f;
}

InstableFloat C()
{  /* C is the lower limit of Y */ 
   double f; 

   f = 0.0; 
   return f;
}

InstableFloat D(InstableFloat X)
{  /* D is the upper limit of Y */
  InstableFloat f, t;
  InstableEnv::open_label(8);
  t = 4 - X * X;
  InstableEnv::close_label();

  f = t;

  return f;
}

InstableFloat ALPHA(InstableFloat X, InstableFloat Y)
{  /* ALPHA is the lower limit for Z */ 
  InstableFloat f,t;
  
  InstableEnv::open_label(9);
  t = X * X + Y * Y;
  InstableEnv::close_label();

  f = t;
 
  return f;
}

InstableFloat BETA()
{  /* BETA is the upper limit for Z */ 
   InstableFloat f; 

   f = 2; 
   return f;
}


void INPUT(double *A, double *B)
{
  //double X; 
  //char AA;
   FILE *INP;
   INP = fopen("triple.txt", "r");
   for(int i = 0; i < 10; i ++){
     fscanf(INP, "%lf %lf", &A[i], &B[i]);
   }
}


