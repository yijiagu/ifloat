
/*
*   CROUT REDUCTION FOR TRIDIAGONAL LINEAR SYSTEMS ALGORITHM 6.7
*
*   To solve the n x n linear system
*
*   E1:  A[1,1] X[1] + A[1,2] X[2]                  = A[1,n+1]
*   E2:  A[2,1] X[1] + A[2,2] X[2] + A[2,3] X[3]    = A[2,n+1]
*   :
*   .
*   E(n):          A[n,n-1] X[n-1] + A[n,n] X[n]    = A[n,n+1]
*
*   INPUT:   the dimension n; the entries of A.
*
*   OUTPUT:  the solution X(1), ..., X(N).
*/


#include <Random.h>

#include<stdio.h>
#include <iostream>
#include <chrono>
using namespace std::chrono;

#define true 1
#define false 0




int main()
{
   float A[10],B[10],C[10],BB[10],Z[10],X[10];
   int N,I,NN,II;

   N = 10;
   int RANDOM_SEED[10] = { 101010,  
			  111010,  110001, 110011,
			   100001, 101001, 101101,110010, 110110,110010};
   
   high_resolution_clock::time_point t_start = high_resolution_clock::now();   

   for(int i = 0; i < 10; i ++){
     Random R = new_Random_seed(RANDOM_SEED[i]);
     double *A_float = RandomVector(N, R);
     double *B_float = RandomVector(N, R);
     double *C_float = RandomVector(N, R);
     double *BB_float = RandomVector(N, R);

     for (I=1; I<=N; I++)
       A[I-1] = A_float[I-1];

     /* the lower sub-diagonal A(I,I-1) is stored
	in B(I), 2 <= I <= n */ 
     for (I=2; I<=N; I++)
       B[I-1] = B_float[I-1];
     
     /* the upper sub-diagonal A(I,I+1) is stored
	in C(I), 1 <= I <= n-1 */ 
     NN = N - 1;
     for (I=1; I<=NN; I++)
       C[I-1] = C_float[I-1];
     
     for (I=1; I<=N; I++)
       BB[I-1] = BB_float[I - 1];
     
      /* Steps 1-3 set up and solve LZ = B   */
      /* STEP 1 */
      /* the entries of U overwrite C and
            the entries of L overwrite A */
      C[0] = C[0] / A[0]; 
      Z[0] = BB[0] / A[0];
      /* STEP 2 */
      for (I=2; I<=NN; I++) {

	A[I-1] = A[I-1] - B[I-1] * C[I-2];

	
	C[I-1] = C[I-1] / A[I-1];
	float t;

	t = BB[I-1]-B[I-1]*Z[I-2];

	Z[I-1] = t/A[I-1];
      }
      /* STEP 3 */

      A[N-1] = A[N-1] - B[N-1] * C[N-2];

      float t;

      t = BB[N-1]-B[N-1]*Z[N-2];

      Z[N-1] = t/A[N-1];
      /* STEP 4 */
      /* STEPS 4, 5 solve UX = Z */
      X[N-1] = Z[N-1];
      /* STEP 5 */
      for (II=1; II<=NN; II++) {
	I = NN - II + 1;

	X[I-1] = Z[I-1] - C[I-1] * X[I];

      }

      Random_delete(R);
      free(A_float);
      free(B_float);
      free(C_float);
      free(BB_float);
   }

   high_resolution_clock::time_point t_end = high_resolution_clock::now();
   auto duration = duration_cast<microseconds>( t_end - t_start).count();
   std::cout << "duration " << duration << std::endl;
   
   return 0;
}
