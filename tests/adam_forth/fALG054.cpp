/*
*   ADAMS-FORTH ORDER PREDICTOR-CORRECTOR ALGORITHM 5.4
*
*   To approximate the solution of the initial value problem
*          y' = f(t,y), a <= t <= b, y(a) = alpha,
*   at N+1 equally spaced points in the interval [a,b].
*
*   INPUT:   endpoints a,b; initial condition alpha; integer N.
*
*   OUTPUT:  approximation w to y at the (N+1) values of t.
*/

#include <stdio.h>
#include <math.h>

#include <iostream>
#include <chrono>
using namespace std::chrono;

#define true 1
#define false 0

float F(float, float);
void RK4(float , float , float ,
	 float *, float *);

void INPUT(int *, double *, double *, double *, int *);

int main(){
   float T[4], W[4];
   double A_float[10], B_float[10], ALPHA_float[10], H;
   int I, N, J, OK;
   //FILE *OUP[1];

   float A, ALPHA, T0, W0;
   
   INPUT(&OK, A_float, B_float, ALPHA_float, &N);
   
   high_resolution_clock::time_point t_start = high_resolution_clock::now();
   for(int i = 0; i < 10; i ++){
     A = A_float[i];
     ALPHA = ALPHA_float[i];

     if (OK) {
       //OUTPUT(OUP);
       /* STEP 1 */
       H = (B_float[i] - A_float[i]) / N;
      
       T[0] = A;
       W[0] = ALPHA;

       /* STEP 2 */
       for (I=1; I<=3; I++) {
         /* STEP 3 AND 4 */
         /* compute starting values using Runge-Kutta method */
         /* given in a subroutine */
         RK4(H, T[I-1], W[I-1], &T[I], &W[I]);

         /* STEP 5 */
       }
       /* STEP 6 */
       for (I=4; I<=N; I++) {
         /* STEP 7 */
         /* T0, W0 will be used in place of t, w resp. */
         T0 = A + I * H;
         /* predict W(I) */
	 float t1;
	 t1 = 55.0*F(T[3],W[3])-59.0*F(T[2],W[2])
	   +37.0*F(T[1],W[1])-9.0*F(T[0],W[0]);
	 
	 
         W0 = W[3] + H*t1/24.0;
         /* correct W(I) */
	 float t2;
	 t2 = 9.0*F(T0,W0)+19.0*F(T[3],W[3])
	   -5.0*F(T[2],W[2])+F(T[1],W[1]);



         W0 = W[3] + H*t2/24.0;
         /* STEP 8 */

         /* STEP 9 */
         /* prepare for next iteration */
         for (J=1; J<=3; J++) {
	   T[J-1] = T[J];
	   W[J-1] = W[J];
         }
         /* STEP 10 */
         T[3] = T0;
         W[3] = W0;
       }
     }
   }

   high_resolution_clock::time_point t_end = high_resolution_clock::now();
   auto duration = duration_cast<microseconds>( t_end - t_start).count();
   std::cout << "duration " << duration << std::endl;
   
   return 0;
}

/*  Change function F for a new problem    */
float F(float T, float Y){
   float f;
   f = Y - T*T + 1.0;   
   return f;
}

void RK4(float H, float T0, float W0,
	 float *T1, float *W1)
{
   float  K1,K2,K3,K4;
   *T1 = T0 + H;
   K1 = H * F(T0, W0);
   
   float t1, t2, t3, t4, t5, t6;
   t1 = T0 + 0.5 * H;
   t2 = W0 + 0.5 * K1;
   K2 = H * F(t1, t2);
   
   t3 = W0 + 0.5 * K2;
   K3 = H * F(t1, t3);
   
   t4 = K2 + K3;
   t6 = W0 + K3;
   K4 = H * F(*T1, t6);
   
   t5 = K1 + 2.0*t4 + K4;
  
   *W1 = W0 + t5 / 6.0;
}

void INPUT(int *OK, double *A, double *B, double *ALPHA, int *N)
{
  //double X; 
  //char AA;

   printf("This is Adams-Bashforth Predictor Corrector Method\n");
   printf("Has the function F been created in the program \n");
   printf("preceding the INPUT procedure?  Enter Y or N.\n");

   *OK = true;
   *N = 5;
   FILE *INP;
   INP = fopen("adam.txt", "r");
   for(int i = 0; i < 10; i ++){
     fscanf(INP, "%lf %lf %lf", &A[i], &B[i], &ALPHA[i]);
   }
}


